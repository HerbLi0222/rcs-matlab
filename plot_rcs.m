% PLOT_RCS  Plot RCS results from a .dat output file
%
%   Reads the simulation output .dat file and auto-detects data type:
%     - 1D cut (ip==1 or it==1):
%         → 2D Cartesian plot (RCS vs angle)
%         → 2D Polar plot (RCS on polar axes)
%     - 2D grid (ip>1 and it>1):
%         → 3D Theta-Phi mesh/surface plot
%         → 3D Spherical polar plot
%
%   Usage:
%     >> plot_rcs
%     (Then select a .dat file from the file dialog)
%
%     >> plot_rcs('results/temp_20260616220000.dat')

function plot_rcs(filePath)
    clc;

    %% Select file if not provided
    if nargin < 1 || isempty(filePath)
        [fileName, pathName] = uigetfile({'*.dat', 'RCS Result Files (*.dat)'}, ...
                                          'Select RCS result file', 'results/');
        if isequal(fileName, 0)
            disp('No file selected. Exiting.');
            return;
        end
        filePath = fullfile(pathName, fileName);
    end

    [~, shortName, ~] = fileparts(filePath);
    fprintf('========================================\n');
    fprintf('  RCS Result Visualization\n');
    fprintf('========================================\n');
    fprintf('  File: %s\n', shortName);

    %% Parse the .dat file
    [theta, Sth, phi, Sph, param, ip, it] = parseRcsFile(filePath);

    fprintf('  Grid: %d phi x %d theta points\n', ip, it);
    fprintf('  Theta: %.1f ~ %.1f deg\n', min(theta(:)), max(theta(:)));
    fprintf('  Phi:   %.1f ~ %.1f deg\n', min(phi(:)), max(phi(:)));

    %% Extract simulation info
    mode = extractParam(param, 'Mode');
    freq = extractParam(param, 'Radar Frequency');
    pol  = extractParam(param, 'Incident wave polarization');
    fprintf('  Mode: %s | Freq: %s GHz | Pol: %s\n', mode, freq, pol);

    %% Auto-detect: 1D cut or 2D grid?
    is1D = (ip == 1) || (it == 1);

    if is1D
        fprintf('\n  >> Detected 1D angular cut\n');
        fprintf('  >> Generating: Cartesian + Polar plots\n\n');
        plot1D(theta, Sth, phi, Sph, ip, it, shortName, mode, freq, pol);
    else
        fprintf('\n  >> Detected 2D theta-phi grid\n');
        fprintf('  >> Generating: Theta-Phi 3D mesh + 3D Spherical polar\n\n');
        plot2D(theta, Sth, phi, Sph, ip, it, shortName, mode, freq, pol);
    end

    fprintf('Done.\n');
end

%% ==================== Parameter Extraction ====================
function val = extractParam(param, key)
    val = '?';
    idx = strfind(param, key);
    if ~isempty(idx)
        rest = param(idx(1)+length(key):end);
        rest = strrep(rest, ':', '');
        rest = strrep(rest, '(', '');
        rest = strrep(rest, ')', '');
        token = strtrim(rest);
        nl = strfind(token, newline);
        if ~isempty(nl)
            token = token(1:nl(1)-1);
        end
        val = strtrim(token);
    end
end

%% ==================== File Parser (.dat format) ====================
function [theta, Sth, phi, Sph, param, ip, it] = parseRcsFile(filePath)
    fid = fopen(filePath, 'r');
    if fid == -1
        error('Cannot open file: %s', filePath);
    end

    content = {};
    while ~feof(fid)
        content{end+1} = fgetl(fid); %#ok<AGROW>
    end
    fclose(fid);

    % Find section markers
    paramStart   = find(contains(content, 'Simulation Parameters:'), 1);
    thetaStart   = find(contains(content, 'Theta (deg):'), 1);
    rcsThetaStart = find(contains(content, 'RCS Theta (dBsm):'), 1);
    phiStart     = find(contains(content, 'Phi (deg):'), 1);
    rcsPhiStart  = find(contains(content, 'RCS Phi (dBsm):'), 1);

    % Extract parameter string
    param = '';
    if ~isempty(paramStart)
        for k = (paramStart+1):(thetaStart-1)
            line = strtrim(content{k});
            if ~isempty(line)
                param = [param line newline]; %#ok<AGROW>
            end
        end
    end

    % Parse each section
    theta = parseMatrixSection(content, thetaStart+1, rcsThetaStart-1);
    Sth   = parseMatrixSection(content, rcsThetaStart+1, phiStart-1);
    phi   = parseMatrixSection(content, phiStart+1, rcsPhiStart-1);
    Sph   = parseMatrixSection(content, rcsPhiStart+1, length(content));

    ip = size(theta, 1);
    it = size(theta, 2);
end

function mat = parseMatrixSection(content, rowStart, rowEnd)
    rows = {};
    for k = rowStart:rowEnd
        if k > length(content), break; end
        line = strtrim(content{k});
        if isempty(line), continue; end
        line = strrep(line, '[', '');
        line = strrep(line, ']', '');
        line = strrep(line, ';', '');
        line = strtrim(line);
        if isempty(line), continue; end
        nums = sscanf(line, '%f');
        if ~isempty(nums)
            rows{end+1} = nums(:)'; %#ok<AGROW>
        end
    end

    if isempty(rows)
        mat = []; return;
    end

    nCols = max(cellfun(@length, rows));
    mat = zeros(length(rows), nCols);
    for r = 1:length(rows)
        rowData = rows{r};
        n = min(length(rowData), nCols);
        mat(r, 1:n) = rowData(1:n);
    end
end

%% ================================================================
%  1D PLOTS: Cartesian + Polar
%% ================================================================
function plot1D(theta, Sth, phi, Sph, ip, it, shortName, mode, freq, pol)

    if ip == 1
        angleDeg  = theta(1, :);
        SthData   = Sth(1, :);
        SphData   = Sph(1, :);
        angleLabel = '\theta (deg)';
        sliceInfo  = sprintf('phi = %.0f°', phi(1,1));
        sweepVar   = 'theta';
    else
        angleDeg  = phi(:, 1)';
        SthData   = Sth(:, 1)';
        SphData   = Sph(:, 1)';
        angleLabel = '\phi (deg)';
        sliceInfo  = sprintf('\\theta = %.0f°', theta(1,1));
        sweepVar   = 'phi';
    end

    % ---- Figure 1: Cartesian ----
    figure('Name', 'RCS - 1D Cartesian', 'NumberTitle', 'off', ...
           'Position', [100 100 750 480]);

    plot(angleDeg, SthData, 'b-', 'LineWidth', 1.8);
    hold on;
    plot(angleDeg, SphData, 'r--', 'LineWidth', 1.8);

    xlabel(angleLabel, 'FontSize', 12);
    ylabel('RCS (dBsm)', 'FontSize', 12);
    title(sprintf('RCS vs %s  [%s]  |  %s  |  %s GHz  |  %s', ...
          sweepVar, sliceInfo, mode, freq, pol), 'FontSize', 11);
    legend('S_{\theta}', 'S_{\phi}', 'Location', 'best', 'FontSize', 11);
    grid on; box on;
    xlim('tight');

    % Annotate peak
    [peakVal, peakIdx] = max(SthData);
    plot(angleDeg(peakIdx), peakVal, 'bv', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    text(angleDeg(peakIdx), peakVal, ...
         sprintf('  Peak: %.1f dBsm @ %.0f°', peakVal, angleDeg(peakIdx)), ...
         'FontSize', 9, 'Color', 'b');

    % ---- Figure 2: Polar ----
    figure('Name', 'RCS - 1D Polar', 'NumberTitle', 'off', ...
           'Position', [200 100 650 550]);

    minAll = min(min(SthData), min(SphData));
    rTh = SthData - minAll + 3;
    rPh = SphData - minAll + 3;

    polarplot(deg2rad(angleDeg), rTh, 'b-', 'LineWidth', 1.8);
    hold on;
    polarplot(deg2rad(angleDeg), rPh, 'r--', 'LineWidth', 1.8);

    title(sprintf('Polar RCS Pattern  [%s]  |  %s  |  %s', ...
          sliceInfo, mode, pol), 'FontSize', 12);
    legend('S_{\theta}', 'S_{\phi}', 'Location', 'best', 'FontSize', 11);

    maxR = max(max(rTh), max(rPh));
    rLevels = round(linspace(0, maxR, 5));
    for lv = rLevels(2:end-1)
        dbVal = lv + minAll - 3;
        text(lv, 0, sprintf('%.0f dBsm', dbVal), 'FontSize', 7, ...
             'HorizontalAlignment', 'center', 'Color', [0.4 0.4 0.4]);
    end
end

%% ================================================================
%  2D PLOTS: Theta-Phi Mesh + 3D Spherical Polar
%% ================================================================
function plot2D(theta, Sth, phi, Sph, ip, it, shortName, mode, freq, pol)

    phiRange   = max(phi(:)) - min(phi(:));
    thetaRange = max(theta(:)) - min(theta(:));
    needPhiWrap   = (phiRange < 350) && (ip > 1);
    needThetaWrap = (thetaRange < 350) && (it > 1);

    % ---- Figure 1: Theta-Phi 3D Mesh ----
    figure('Name', 'RCS - 3D Mesh', 'NumberTitle', 'off', ...
           'Position', [50 80 1100 480]);

    subplot(1, 2, 1);
    plotThetaPhiMesh(theta, phi, Sth, ip, it, needPhiWrap, needThetaWrap);
    title(sprintf('S_{\\theta} (dBsm)  |  %s  |  %s GHz', mode, freq), 'FontSize', 12);
    xlabel('\theta (deg)'); ylabel('\phi (deg)'); zlabel('RCS (dBsm)');
    colorbar; grid on;

    subplot(1, 2, 2);
    plotThetaPhiMesh(theta, phi, Sph, ip, it, needPhiWrap, needThetaWrap);
    title(sprintf('S_{\\phi} (dBsm)  |  %s  |  %s', mode, pol), 'FontSize', 12);
    xlabel('\theta (deg)'); ylabel('\phi (deg)'); zlabel('RCS (dBsm)');
    colorbar; grid on;

    sgtitle(sprintf('RCS Theta-Phi Pattern  —  %s', shortName), 'FontSize', 13);

    % ---- Figure 2: 3D Spherical Polar ----
    figure('Name', 'RCS - 3D Spherical', 'NumberTitle', 'off', ...
           'Position', [50 50 1100 500]);

    minBoth = min(min(Sth(:)), min(Sph(:)));

    subplot(1, 2, 1);
    plotSphericalPolar(theta, phi, Sth, minBoth, ip, it);
    title(sprintf('S_{\\theta} (dBsm)  |  %s  |  %s GHz', mode, freq), 'FontSize', 12);
    colorbar;

    subplot(1, 2, 2);
    plotSphericalPolar(theta, phi, Sph, minBoth, ip, it);
    title(sprintf('S_{\\phi} (dBsm)  |  %s  |  %s', mode, pol), 'FontSize', 12);
    colorbar;

    sgtitle(sprintf('3D Spherical RCS Pattern  —  %s', shortName), 'FontSize', 13);
end

%% --- Theta-Phi mesh helper ---
function plotThetaPhiMesh(theta, phi, data, ip, it, needPhiWrap, needThetaWrap)
    phiDeg   = phi;
    thetaDeg = theta;
    plotData = data;

    if needPhiWrap && ip > 1
        phiDeg   = [phiDeg; phiDeg(1,:) + 360];
        plotData = [plotData; plotData(1,:)];
    end
    if needThetaWrap && it > 1
        thetaDeg = [thetaDeg, thetaDeg(:,1)];
        plotData = [plotData, plotData(:,1)];
    end

    [T, P] = meshgrid(thetaDeg(1,:), phiDeg(:,1));
    if size(plotData,1) ~= size(P,1)
        plotData = plotData';
    end

    surf(T, P, plotData, 'EdgeColor', 'interp', 'FaceColor', 'interp');
    colormap(jet);
    view(135, 35);
end

%% --- 3D Spherical polar helper ---
function plotSphericalPolar(theta, phi, RCS, minRCS, ip, it)
    thVec = deg2rad(theta(1, :));
    phVec = deg2rad(phi(:, 1));

    phFull  = [phVec; phVec(1) + 2*pi];
    rcsData = RCS;
    rcsData(end+1, :) = RCS(1, :);

    [ThetaGrid, PhiGrid] = meshgrid(thVec, phFull);

    rMin  = 0.3;
    rData = RCS - minRCS + rMin;
    rData(end+1, :) = rData(1, :);

    X = rData .* sin(ThetaGrid) .* cos(PhiGrid);
    Y = rData .* sin(ThetaGrid) .* sin(PhiGrid);
    Z = rData .* cos(ThetaGrid);

    surf(X, Y, Z, rcsData, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
    hold on;

    [Sx, Sy, Sz] = sphere(40);
    Sx = Sx * rMin; Sy = Sy * rMin; Sz = Sz * rMin;
    mesh(Sx, Sy, Sz, 'FaceColor', 'none', 'EdgeColor', [0.6 0.6 0.6], ...
         'EdgeAlpha', 0.25, 'LineWidth', 0.3);

    mx = max(abs([X(:); Y(:); Z(:)])) * 1.2;
    line([-mx mx], [0 0], [0 0], 'Color', 'k', 'LineWidth', 0.5);
    line([0 0], [-mx mx], [0 0], 'Color', 'k', 'LineWidth', 0.5);
    line([0 0], [0 0], [-mx mx], 'Color', 'k', 'LineWidth', 0.5);
    text(mx*1.05, 0, 0, 'X', 'FontSize', 10);
    text(0, mx*1.05, 0, 'Y', 'FontSize', 10);
    text(0, 0, mx*1.05, 'Z', 'FontSize', 10);

    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal; axis tight;
    view(135, 25);
    lighting gouraud;
    camlight('headlight');
    colormap(jet);
end
