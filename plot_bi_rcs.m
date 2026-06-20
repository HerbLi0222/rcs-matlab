% PLOT_BI_RCS  Plot bistatic RCS 2D images from a .dat output file
%
%   Reads the bistatic simulation output .dat file and generates
%   2D visualization of the bistatic scattering pattern:
%     - 2D Theta-Phi color mesh (RCS vs observation angles)
%     - U-V contour plot (direction cosine space)
%     - Polar heatmap of the scattering hemisphere
%
%   Bistatic data: fixed incident direction, RCS over a 2D observation
%   grid (theta, phi). The file format is the same as monostatic .dat.
%
%   Usage:
%     >> plot_bi_rcs
%     (Then select a .dat file from the file dialog)
%
%     >> plot_bi_rcs('results/temp_20260619222000.dat')
%
%   See also: plot_mono_rcs, main_bistatic

function plot_bi_rcs(filePath)
    clc;

    %% Select file if not provided
    if nargin < 1 || isempty(filePath)
        [fileName, pathName] = uigetfile({'*.dat', 'RCS Result Files (*.dat)'}, ...
                                          'Select bistatic RCS result file', 'results/');
        if isequal(fileName, 0)
            disp('No file selected. Exiting.');
            return;
        end
        filePath = fullfile(pathName, fileName);
    end

    [~, shortName, ~] = fileparts(filePath);
    fprintf('========================================\n');
    fprintf('  Bistatic RCS 2D Visualization\n');
    fprintf('========================================\n');
    fprintf('  File: %s\n', shortName);

    %% Parse the .dat file
    [theta, Sth, phi, Sph, param, ip, it] = parseBiRcsFile(filePath);

    fprintf('  Grid: %d phi x %d theta points\n', ip, it);
    fprintf('  Theta: %.1f ~ %.1f deg\n', min(theta(:)), max(theta(:)));
    fprintf('  Phi:   %.1f ~ %.1f deg\n', min(phi(:)), max(phi(:)));

    %% Extract simulation info
    mode = extractBiParam(param, 'Mode');
    freq = extractBiParam(param, 'Radar Frequency');
    pol  = extractBiParam(param, 'Incident wave polarization');
    fprintf('  Mode: %s | Freq: %s GHz | Pol: %s\n', mode, freq, pol);

    %% Validate: bistatic data should be 2D
    is1D = (ip == 1) || (it == 1);
    if is1D
        fprintf('\n  >> Warning: Data appears to be a 1D cut, not full 2D bistatic.\n');
        fprintf('  >> Use plot_mono_rcs for 1D cuts.\n');
        fprintf('  >> Attempting 1D plot anyway...\n\n');
        plotBi1D(theta, Sth, phi, Sph, ip, it, shortName, mode, freq, pol);
        return;
    end

    fprintf('\n  >> Detected 2D bistatic pattern\n');

    %% Extract monostatic point from bistatic data
    % Monostatic condition: observation direction = incidence direction
    monoInfo = extractMonoPoint(theta, phi, Sth, Sph, param);

    fprintf('\n  >> Monostatic point (obs = inc):\n');
    fprintf('     \\theta = %.1f°, \\phi = %.1f°  →  S_\\theta = %.2f dBsm, S_\\phi = %.2f dBsm\n', ...
            monoInfo.thMono, monoInfo.phMono, monoInfo.SthMono, monoInfo.SphMono);

    fprintf('  >> Generating: Theta-Phi mesh + U-V contour + Polar heatmap + Monostatic\n\n');

    %% ---- Figure 1: Theta-Phi 2D Color Mesh ----
    plotBiThetaPhi(theta, phi, Sth, Sph, ip, it, shortName, mode, freq, pol, param, monoInfo);

    %% ---- Figure 2: U-V Direction Cosine Contour ----
    plotBiUV(theta, phi, Sth, Sph, ip, it, shortName, mode, freq, pol, param, monoInfo);

    %% ---- Figure 3: Polar Heatmap ----
    plotBiPolarHeatmap(theta, phi, Sth, ip, it, shortName, mode, freq, pol, param, monoInfo);

    %% ---- Figure 4: Monostatic RCS extracted from bistatic ----
    plotBiMono(theta, phi, Sth, Sph, ip, it, shortName, mode, freq, pol, param, monoInfo);

    fprintf('Done.\n');
end

%% ==================== Parameter Extraction ====================
function val = extractBiParam(param, key)
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
function [theta, Sth, phi, Sph, param, ip, it] = parseBiRcsFile(filePath)
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
    paramStart    = find(contains(content, 'Simulation Parameters:'), 1);
    thetaStart    = find(contains(content, 'Theta (deg):'), 1);
    rcsThetaStart = find(contains(content, 'RCS Theta (dBsm):'), 1);
    phiStart      = find(contains(content, 'Phi (deg):'), 1);
    rcsPhiStart   = find(contains(content, 'RCS Phi (dBsm):'), 1);

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
    theta = parseBiMatrixSection(content, thetaStart+1, rcsThetaStart-1);
    Sth   = parseBiMatrixSection(content, rcsThetaStart+1, phiStart-1);
    phi   = parseBiMatrixSection(content, phiStart+1, rcsPhiStart-1);
    Sph   = parseBiMatrixSection(content, rcsPhiStart+1, length(content));

    ip = size(theta, 1);
    it = size(theta, 2);
end

function mat = parseBiMatrixSection(content, rowStart, rowEnd)
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
%  1D fallback plot (shouldn't normally be used for bistatic)
%% ================================================================
function plotBi1D(theta, Sth, phi, Sph, ip, it, shortName, mode, freq, pol)

    if ip == 1
        angleDeg  = theta(1, :);
        SthData   = Sth(1, :);
        SphData   = Sph(1, :);
        angleLabel = '\theta (deg)';
        sliceInfo  = sprintf('\\phi = %.0f°', phi(1,1));
    else
        angleDeg  = phi(:, 1)';
        SthData   = Sth(:, 1)';
        SphData   = Sph(:, 1)';
        angleLabel = '\phi (deg)';
        sliceInfo  = sprintf('\\theta = %.0f°', theta(1,1));
    end

    figure('Name', 'Bistatic RCS - 1D Cut', 'NumberTitle', 'off', ...
           'Position', [100 100 750 480]);

    plot(angleDeg, SthData, 'b-', 'LineWidth', 1.8);
    hold on;
    plot(angleDeg, SphData, 'r--', 'LineWidth', 1.8);

    xlabel(angleLabel, 'FontSize', 12);
    ylabel('RCS (dBsm)', 'FontSize', 12);
    title(sprintf('Bistatic RCS  [%s]  |  %s  |  %s GHz  |  %s', ...
          sliceInfo, mode, freq, pol), 'FontSize', 11);
    legend('S_{\theta}', 'S_{\phi}', 'Location', 'best', 'FontSize', 11);
    grid on; box on;
    xlim('tight');
end

%% ================================================================
%  Figure 1: Theta-Phi 2D Color Mesh
%    Full bistatic scattering pattern over observation angles
%% ================================================================
function plotBiThetaPhi(theta, phi, Sth, Sph, ip, it, shortName, mode, freq, pol, param, monoInfo)
    % Extract incidence info for title
    incTheta = extractBiParam(param, 'Start theta');
    incPhi   = extractBiParam(param, 'Start phi');

    figure('Name', 'Bistatic RCS - Theta-Phi Pattern', 'NumberTitle', 'off', ...
           'Position', [50 80 1200 500]);

    % ---- S_theta subplot ----
    subplot(1, 2, 1);
    plotBiMesh(theta, phi, Sth, ip, it);
    hold on;
    % Mark monostatic point
    plot3(monoInfo.thMono, monoInfo.phMono, monoInfo.SthMono, ...
          'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'k', ...
          'DisplayName', sprintf('Mono: %.1f dBsm', monoInfo.SthMono));
    hold off;
    title(sprintf('S_{\\theta} (dBsm)  |  %s  |  %s GHz', mode, freq), 'FontSize', 12);
    xlabel('\theta_{obs} (deg)', 'FontSize', 11);
    ylabel('\phi_{obs} (deg)', 'FontSize', 11);
    zlabel('RCS (dBsm)', 'FontSize', 11);
    cbar = colorbar; ylabel(cbar, 'dBsm', 'FontSize', 10);
    grid on; legend('Location', 'best', 'FontSize', 8);

    % ---- S_phi subplot ----
    subplot(1, 2, 2);
    plotBiMesh(theta, phi, Sph, ip, it);
    hold on;
    plot3(monoInfo.thMono, monoInfo.phMono, monoInfo.SphMono, ...
          'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'k', ...
          'DisplayName', sprintf('Mono: %.1f dBsm', monoInfo.SphMono));
    hold off;
    title(sprintf('S_{\\phi} (dBsm)  |  %s  |  %s', mode, pol), 'FontSize', 12);
    xlabel('\theta_{obs} (deg)', 'FontSize', 11);
    ylabel('\phi_{obs} (deg)', 'FontSize', 11);
    zlabel('RCS (dBsm)', 'FontSize', 11);
    cbar = colorbar; ylabel(cbar, 'dBsm', 'FontSize', 10);
    grid on; legend('Location', 'best', 'FontSize', 8);

    sgtitle(sprintf('Bistatic RCS Pattern  —  %s  |  Incidence: \\theta_i=%s°, \\phi_i=%s°  |  ● = Monostatic', ...
            shortName, incTheta, incPhi), 'FontSize', 13);
end

function plotBiMesh(theta, phi, data, ip, it)
    % Get 1D angle vectors
    thVec = theta(1, :);   % 1 x it
    phVec = phi(:, 1);     % ip x 1

    phiRange = max(phVec) - min(phVec);
    needPhiWrap = (phiRange < 350) && (ip > 1);

    if needPhiWrap
        phVec = [phVec; phVec(1) + 360];
        plotData = [data; data(1, :)];
    else
        plotData = data;
    end

    [T, P] = meshgrid(thVec, phVec);

    % Ensure data orientation matches mesh
    if size(plotData, 1) ~= size(P, 1)
        plotData = plotData';
    end

    surf(T, P, plotData, 'EdgeColor', 'interp', 'FaceColor', 'interp');
    colormap(jet);
    view(135, 35);
end

%% ================================================================
%  Figure 2: U-V Direction Cosine Contour Plot
%    Scattering pattern in direction cosine (U,V) space
%    U = sin(theta)*cos(phi), V = sin(theta)*sin(phi)
%% ================================================================
function plotBiUV(theta, phi, Sth, Sph, ip, it, shortName, mode, freq, pol, param, monoInfo)
    incTheta = extractBiParam(param, 'Start theta');
    incPhi   = extractBiParam(param, 'Start phi');

    % Convert to U-V space
    thRad = deg2rad(theta);
    phRad = deg2rad(phi);
    U = sin(thRad) .* cos(phRad);
    V = sin(thRad) .* sin(phRad);

    % Mark incidence direction (also the monostatic direction)
    thIncRad = deg2rad(str2double(incTheta));
    phIncRad = deg2rad(str2double(incPhi));
    U_inc = sin(thIncRad) * cos(phIncRad);
    V_inc = sin(thIncRad) * sin(phIncRad);

    % Monostatic point in U-V space
    thMonoRad = deg2rad(monoInfo.thMono);
    phMonoRad = deg2rad(monoInfo.phMono);
    U_mono = sin(thMonoRad) * cos(phMonoRad);
    V_mono = sin(thMonoRad) * sin(phMonoRad);

    figure('Name', 'Bistatic RCS - UV Contour', 'NumberTitle', 'off', ...
           'Position', [50 50 1200 500]);

    % ---- S_theta ----
    subplot(1, 2, 1);
    contourf(U, V, Sth, 40, 'LineColor', 'none');
    hold on;
    plot(U_inc, V_inc, 'wx', 'MarkerSize', 14, 'LineWidth', 2.5, ...
         'DisplayName', sprintf('Inc = Mono: \\theta=%.0f°,\\phi=%.0f°', ...
                                str2double(incTheta), str2double(incPhi)));
    plot(U_mono, V_mono, 'ko', 'MarkerSize', 12, 'LineWidth', 2, ...
         'DisplayName', sprintf('Mono S_\\theta=%.1f dBsm', monoInfo.SthMono));
    colormap(jet);
    cbar = colorbar; ylabel(cbar, 'dBsm', 'FontSize', 10);
    xlabel('U = sin\theta cos\phi', 'FontSize', 11);
    ylabel('V = sin\theta sin\phi', 'FontSize', 11);
    title(sprintf('S_{\\theta}  |  %s  |  %s GHz', mode, freq), 'FontSize', 12);
    axis equal; axis([-1 1 -1 1]);
    grid on; box on;
    legend('Location', 'best', 'FontSize', 9);

    % ---- S_phi ----
    subplot(1, 2, 2);
    contourf(U, V, Sph, 40, 'LineColor', 'none');
    hold on;
    plot(U_inc, V_inc, 'wx', 'MarkerSize', 14, 'LineWidth', 2.5, ...
         'DisplayName', sprintf('Inc = Mono: \\theta=%.0f°,\\phi=%.0f°', ...
                                str2double(incTheta), str2double(incPhi)));
    plot(U_mono, V_mono, 'ko', 'MarkerSize', 12, 'LineWidth', 2, ...
         'DisplayName', sprintf('Mono S_\\phi=%.1f dBsm', monoInfo.SphMono));
    colormap(jet);
    cbar = colorbar; ylabel(cbar, 'dBsm', 'FontSize', 10);
    xlabel('U = sin\theta cos\phi', 'FontSize', 11);
    ylabel('V = sin\theta sin\phi', 'FontSize', 11);
    title(sprintf('S_{\\phi}  |  %s  |  %s', mode, pol), 'FontSize', 12);
    axis equal; axis([-1 1 -1 1]);
    grid on; box on;
    legend('Location', 'best', 'FontSize', 9);

    sgtitle(sprintf('Bistatic UV Pattern  —  %s  |  X = Incidence  |  O = Monostatic', ...
            shortName), 'FontSize', 13);
end

%% ================================================================
%  Figure 3: Polar Heatmap (Scattering Hemisphere)
%    Observation hemisphere with RCS mapped on polar axes
%% ================================================================
function plotBiPolarHeatmap(theta, phi, Sth, ip, it, shortName, mode, freq, pol, param, monoInfo)
    incTheta = extractBiParam(param, 'Start theta');
    incPhi   = extractBiParam(param, 'Start phi');

    % Normalize data for colormap
    SthNorm = Sth;
    Smax = max(SthNorm(:));
    Smin = min(SthNorm(:));

    % Build polar grid
    thVec = unique(theta(1, :));  % theta values
    phVec = unique(phi(:, 1));    % phi values

    % Create a higher-resolution grid for smooth polar display
    nTh = min(length(thVec), 181);
    nPh = min(length(phVec), 361);
    thHi = linspace(min(thVec), max(thVec), nTh);
    phHi = linspace(min(phVec), max(phVec), nPh);

    [TH, PH] = meshgrid(thHi, phHi);
    [TH_orig, PH_orig] = meshgrid(thVec, phVec);

    % Interpolate to high-res grid
    if size(Sth, 1) == length(phVec) && size(Sth, 2) == length(thVec)
        SthInterp = interp2(TH_orig, PH_orig, Sth, TH, PH, 'linear');
    else
        SthInterp = interp2(TH_orig, PH_orig, Sth', TH, PH, 'linear');
    end

    figure('Name', 'Bistatic RCS - Polar Heatmap', 'NumberTitle', 'off', ...
           'Position', [100 60 700 600]);

    % Use pcolor on polar-transformed coordinates
    rho = TH;  % theta as radial coordinate (0 to 180 deg)
    [Rho, PhiGrid] = meshgrid(thHi, deg2rad(phHi));

    % Map to Cartesian for visualization
    XX = Rho .* cos(PhiGrid);
    YY = Rho .* sin(PhiGrid);

    pcolor(XX, YY, SthInterp);
    shading interp;
    colormap(jet);
    cbar = colorbar; ylabel(cbar, 'dBsm', 'FontSize', 10);
    caxis([Smin, Smax]);

    hold on;

    % Mark incidence direction (also monostatic direction)
    thInc = str2double(incTheta);
    phInc = str2double(incPhi);
    xi = thInc * cos(deg2rad(phInc));
    yi = thInc * sin(deg2rad(phInc));
    plot(xi, yi, 'wx', 'MarkerSize', 14, 'LineWidth', 2.5);

    % Mark monostatic point (same location, but with circle)
    xm = monoInfo.thMono * cos(deg2rad(monoInfo.phMono));
    ym = monoInfo.thMono * sin(deg2rad(monoInfo.phMono));
    plot(xm, ym, 'ko', 'MarkerSize', 12, 'LineWidth', 1.8);

    % Draw rings at 30° intervals
    for r = 30:30:180
        tRing = linspace(0, 2*pi, 361);
        plot(r * cos(tRing), r * sin(tRing), 'w-', 'LineWidth', 0.3, ...
             'Color', [1 1 1 0.3]);
        if r <= max(thVec)
            text(r, 0, sprintf('%d°', r), 'FontSize', 7, 'Color', 'w', ...
                 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
        end
    end

    % Draw radial spokes
    for a = 0:45:315
        aRad = deg2rad(a);
        plot([0 180*cos(aRad)], [0 180*sin(aRad)], 'w-', 'LineWidth', 0.3, ...
             'Color', [1 1 1 0.3]);
        text(185*cos(aRad), 185*sin(aRad), sprintf('%d°', a), 'FontSize', 8, ...
             'Color', 'w', 'HorizontalAlignment', 'center');
    end

    xlabel('X = \theta cos\phi (°)', 'FontSize', 11);
    ylabel('Y = \theta sin\phi (°)', 'FontSize', 11);
    title(sprintf('Bistatic S_{\\theta} Polar Map  —  %s  |  %s GHz  |  %s\nInc = Mono: \\theta=%s°, \\phi=%s° (X) | O = Mono: %.1f dBsm', ...
          shortName, freq, pol, incTheta, incPhi, monoInfo.SthMono), 'FontSize', 12);
    axis equal; axis tight;
    set(gca, 'Color', [0.15 0.15 0.15]);
end

%% ================================================================
%  Monostatic Point Extraction from Bistatic Data
%    Finds the observation angle closest to the incidence direction.
%    In bistatic mode, monostatic RCS = RCS(θ_obs = θ_inc, φ_obs = φ_inc).
%% ================================================================
function monoInfo = extractMonoPoint(theta, phi, Sth, Sph, param)
    % Get incidence direction from parameters
    incTheta = str2double(extractBiParam(param, 'Start theta'));
    incPhi   = str2double(extractBiParam(param, 'Start phi'));

    % Get observation angle vectors
    thVec = theta(1, :);   % theta observation angles
    phVec = phi(:, 1);     % phi observation angles

    % Find closest observation angles to incidence direction
    [~, thIdx] = min(abs(thVec - incTheta));
    [~, phIdx] = min(abs(phVec - incPhi));

    thMono = thVec(thIdx);
    phMono = phVec(phIdx);

    % Verify data orientation and extract RCS
    if size(Sth, 1) == length(phVec) && size(Sth, 2) == length(thVec)
        SthMono = Sth(phIdx, thIdx);
        SphMono = Sph(phIdx, thIdx);
    else
        SthMono = Sth(thIdx, phIdx);
        SphMono = Sph(thIdx, phIdx);
    end

    % Pack results
    monoInfo = struct();
    monoInfo.thMono   = thMono;
    monoInfo.phMono   = phMono;
    monoInfo.thIdx    = thIdx;
    monoInfo.phIdx    = phIdx;
    monoInfo.SthMono  = SthMono;
    monoInfo.SphMono  = SphMono;
    monoInfo.incTheta = incTheta;
    monoInfo.incPhi   = incPhi;
end

%% ================================================================
%  Figure 4: Monostatic RCS Extraction from Bistatic Data
%    Shows:
%      - 1D cuts through the monostatic point (θ-cut and φ-cut)
%      - The monostatic point highlighted on each cut
%      - Summary panel with extracted monostatic RCS values
%% ================================================================
function plotBiMono(theta, phi, Sth, Sph, ip, it, shortName, mode, freq, pol, param, monoInfo)
    incTheta = extractBiParam(param, 'Start theta');
    incPhi   = extractBiParam(param, 'Start phi');

    thVec = theta(1, :);
    phVec = phi(:, 1);

    %% Extract cut data
    % Theta-cut at phi = phi_mono
    if size(Sth, 1) == length(phVec) && size(Sth, 2) == length(thVec)
        SthThCut = Sth(monoInfo.phIdx, :);
        SphThCut = Sph(monoInfo.phIdx, :);
        SthPhCut = Sth(:, monoInfo.thIdx);
        SphPhCut = Sph(:, monoInfo.thIdx);
    else
        SthThCut = Sth(:, monoInfo.thIdx)';
        SphThCut = Sph(:, monoInfo.thIdx)';
        SthPhCut = Sth(monoInfo.phIdx, :);
        SphPhCut = Sph(monoInfo.phIdx, :);
    end

    figure('Name', 'Bistatic → Monostatic Extraction', 'NumberTitle', 'off', ...
           'Position', [50 80 1100 750]);

    %% -- Subplot 1: Theta-cut at phi = phi_mono --
    subplot(2, 3, 1);
    plot(thVec, SthThCut, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(thVec, SphThCut, 'r--', 'LineWidth', 1.5);
    xline(monoInfo.thMono, 'k--', 'LineWidth', 1.2);
    plot(monoInfo.thMono, monoInfo.SthMono, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    plot(monoInfo.thMono, monoInfo.SphMono, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    hold off;
    xlabel('\theta_{obs} (deg)', 'FontSize', 10);
    ylabel('RCS (dBsm)', 'FontSize', 10);
    title(sprintf('\\theta-Cut at \\phi = %.1f°', monoInfo.phMono), 'FontSize', 11);
    legend('S_{\theta}', 'S_{\phi}', 'Location', 'best', 'FontSize', 8);
    grid on; box on;

    %% -- Subplot 2: Phi-cut at theta = theta_mono --
    subplot(2, 3, 2);
    plot(phVec, SthPhCut, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(phVec, SphPhCut, 'r--', 'LineWidth', 1.5);
    xline(monoInfo.phMono, 'k--', 'LineWidth', 1.2);
    plot(monoInfo.phMono, monoInfo.SthMono, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    plot(monoInfo.phMono, monoInfo.SphMono, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    hold off;
    xlabel('\phi_{obs} (deg)', 'FontSize', 10);
    ylabel('RCS (dBsm)', 'FontSize', 10);
    title(sprintf('\\phi-Cut at \\theta = %.1f°', monoInfo.thMono), 'FontSize', 11);
    legend('S_{\theta}', 'S_{\phi}', 'Location', 'best', 'FontSize', 8);
    grid on; box on;

    %% -- Subplot 3: Summary panel --
    subplot(2, 3, 3);
    axis off;
    text(0.05, 0.95, 'MONOSTATIC RCS', 'FontSize', 14, 'FontWeight', 'bold', ...
         'Units', 'normalized', 'Color', [0.2 0.2 0.6]);
    text(0.05, 0.82, 'Extracted from Bistatic Data', 'FontSize', 10, ...
         'Units', 'normalized', 'Color', [0.4 0.4 0.4]);

    yPos = 0.65;
    dy = 0.07;
    text(0.05, yPos, sprintf('File: %s', shortName), 'FontSize', 9, ...
         'Units', 'normalized', 'Interpreter', 'none'); yPos = yPos - dy;
    text(0.05, yPos, sprintf('Mode: %s | Freq: %s GHz', mode, freq), 'FontSize', 9, ...
         'Units', 'normalized'); yPos = yPos - dy;
    text(0.05, yPos, sprintf('Polarization: %s', pol), 'FontSize', 9, ...
         'Units', 'normalized'); yPos = yPos - dy - 0.02;

    text(0.05, yPos, 'Incidence (= Monostatic) Direction:', 'FontSize', 9, ...
         'Units', 'normalized', 'FontWeight', 'bold'); yPos = yPos - dy;
    text(0.08, yPos, sprintf('\\theta_i = %.1f°   \\phi_i = %.1f°', ...
         str2double(incTheta), str2double(incPhi)), 'FontSize', 10, ...
         'Units', 'normalized', 'Color', [0.2 0.2 0.6]); yPos = yPos - dy - 0.02;

    text(0.05, yPos, 'Observed Monostatic RCS:', 'FontSize', 9, ...
         'Units', 'normalized', 'FontWeight', 'bold'); yPos = yPos - dy;
    text(0.08, yPos, sprintf('S_{\\theta} = %.2f dBsm', monoInfo.SthMono), ...
         'FontSize', 12, 'Units', 'normalized', 'Color', 'b', 'FontWeight', 'bold');
    yPos = yPos - dy;
    text(0.08, yPos, sprintf('S_{\\phi} = %.2f dBsm', monoInfo.SphMono), ...
         'FontSize', 12, 'Units', 'normalized', 'Color', 'r', 'FontWeight', 'bold');
    yPos = yPos - dy - 0.02;

    % Linear value
    SthLin = 10^(monoInfo.SthMono / 10);
    SphLin = 10^(monoInfo.SphMono / 10);
    text(0.08, yPos, sprintf('(Linear: S_{\\theta} = %.4f m^2, S_{\\phi} = %.4f m^2)', ...
         SthLin, SphLin), 'FontSize', 8, 'Units', 'normalized', ...
         'Color', [0.4 0.4 0.4]);

    %% -- Subplot 4: Zoomed theta-cut around monostatic point --
    subplot(2, 3, 4);
    % Zoom to ±20° around monostatic theta
    zoomHalf = 20;
    zoomMask = abs(thVec - monoInfo.thMono) <= zoomHalf;
    plot(thVec(zoomMask), SthThCut(zoomMask), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(thVec(zoomMask), SphThCut(zoomMask), 'r--', 'LineWidth', 1.5);
    xline(monoInfo.thMono, 'k-', 'LineWidth', 1.5);
    plot(monoInfo.thMono, monoInfo.SthMono, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    plot(monoInfo.thMono, monoInfo.SphMono, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    hold off;
    xlabel('\theta_{obs} (deg)', 'FontSize', 10);
    ylabel('RCS (dBsm)', 'FontSize', 10);
    title(sprintf('Zoom \\theta-Cut ±%d° at \\phi=%.1f°', zoomHalf, monoInfo.phMono), 'FontSize', 11);
    legend('S_{\theta}', 'S_{\phi}', 'Location', 'best', 'FontSize', 8);
    grid on; box on;

    %% -- Subplot 5: Zoomed phi-cut around monostatic point --
    subplot(2, 3, 5);
    zoomMask = abs(phVec - monoInfo.phMono) <= zoomHalf;
    plot(phVec(zoomMask), SthPhCut(zoomMask), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(phVec(zoomMask), SphPhCut(zoomMask), 'r--', 'LineWidth', 1.5);
    xline(monoInfo.phMono, 'k-', 'LineWidth', 1.5);
    plot(monoInfo.phMono, monoInfo.SthMono, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    plot(monoInfo.phMono, monoInfo.SphMono, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    hold off;
    xlabel('\phi_{obs} (deg)', 'FontSize', 10);
    ylabel('RCS (dBsm)', 'FontSize', 10);
    title(sprintf('Zoom \\phi-Cut ±%d° at \\theta=%.1f°', zoomHalf, monoInfo.thMono), 'FontSize', 11);
    legend('S_{\theta}', 'S_{\phi}', 'Location', 'best', 'FontSize', 8);
    grid on; box on;

    %% -- Subplot 6: Comparison bar chart --
    subplot(2, 3, 6);
    barData = [monoInfo.SthMono, monoInfo.SphMono];
    b = bar([1 2], barData, 0.5);
    b.FaceColor = 'flat';
    b.CData(1, :) = [0.2 0.4 0.8];
    b.CData(2, :) = [0.8 0.3 0.3];
    set(gca, 'XTickLabel', {'S_{\theta}', 'S_{\phi}'}, 'FontSize', 11);
    ylabel('RCS (dBsm)', 'FontSize', 10);
    title(sprintf('Monostatic RCS  @  \\theta=%.1f°, \\phi=%.1f°', ...
          monoInfo.thMono, monoInfo.phMono), 'FontSize', 11);
    grid on; box on;
    % Add value labels
    for i = 1:2
        text(i, barData(i) + 0.5, sprintf('%.2f', barData(i)), ...
             'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
    end

    sgtitle(sprintf('Monostatic RCS Extracted from Bistatic Data  —  %s', ...
            shortName), 'FontSize', 13, 'Interpreter', 'none');
end
