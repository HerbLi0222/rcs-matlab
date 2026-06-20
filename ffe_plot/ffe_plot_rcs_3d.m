% FFE_PLOT_RCS_3D  3D RCS pattern visualization (远场三维RCS方向图可视化)
%
%   Creates 3D visualizations of 2D angular-grid RCS data:
%     - 3D Spherical surface plot (RCS mapped onto unit sphere)
%     - Theta-Phi 3D mesh/surface (waterfall-style)
%     - 2D Heatmap / pcolor with contour overlay
%
%   For 1D cut data, the spherical and mesh plots are skipped and only
%   the heatmap-style 1D profile is shown.
%
%   Input:
%       data    - data struct from ffe_load_data (must be 2D grid for full 3D)
%       options - (optional) struct with fields:
%           .component     - 'theta' (default), 'phi', or 'both'
%           .plotType      - 'all' (default), 'spherical', 'mesh', 'heatmap'
%           .colormap      - colormap name: 'jet' (default), 'parula', 'hot', 'turbo'
%           .viewAngle     - [az, el] for 3D view (default = [135, 30])
%           .lighting      - true/false for gouraud lighting (default = true)
%           .showSTL       - true/false to overlay STL wireframe (default = false)
%           .stlModel      - path to STL model file (if showSTL is true)
%           .saveFigs      - true/false to save figures (default = false)
%           .saveDir       - directory for saved figures
%           .figPrefix     - prefix for saved figure filenames
%
%   Usage:
%     >> data = ffe_load_data('results/temp_20260616220000.dat');
%     >> ffe_plot_rcs_3d(data);
%     >> ffe_plot_rcs_3d(data, struct('plotType', 'spherical', 'viewAngle', [45, 25]));
%
%   See also: ffe_main, ffe_load_data, ffe_plot_rcs_2d, ffe_utils

function ffe_plot_rcs_3d(data, options)

    %% ---- Default options ----
    if nargin < 2, options = struct(); end
    opt = setDefaults(options, data);

    %% ---- Validate ----
    if ~data.is2D
        warning('ffe_plot_rcs_3d:Not2D', ...
                'Data is a 1D cut. Only heatmap (1D profile) will be plotted.');
    end

    %% ---- Prepare data for selected component(s) ----
    components = {};
    switch lower(opt.component)
        case 'theta'
            components = {'theta'};
        case 'phi'
            components = {'phi'};
        case 'both'
            components = {'theta', 'phi'};
    end

    for cIdx = 1:length(components)
        comp = components{cIdx};
        switch comp
            case 'theta'
                rcsData = data.Sth;
                compLabel = 'S_{\theta}';
            case 'phi'
                rcsData = data.Sph;
                compLabel = 'S_{\phi}';
        end

        %% ---- Spherical 3D plot ----
        if ismember(opt.plotType, {'all', 'spherical'}) && data.is2D
            plotSpherical(data, rcsData, compLabel, cIdx, length(components), opt);
        end

        %% ---- Theta-Phi mesh plot ----
        if ismember(opt.plotType, {'all', 'mesh'}) && data.is2D
            plotThetaPhiMesh(data, rcsData, compLabel, cIdx, length(components), opt);
        end

        %% ---- 2D heatmap / contour ----
        if ismember(opt.plotType, {'all', 'heatmap'})
            plotHeatmap(data, rcsData, compLabel, cIdx, length(components), opt);
        end
    end
end

%% ========================================================================
%  SETDEFAULTS  Fill in default option values
% ========================================================================
function opt = setDefaults(options, data)
    opt.component = getOpt(options, 'component', 'theta');
    opt.plotType  = getOpt(options, 'plotType', 'all');
    opt.colormap  = getOpt(options, 'colormap', 'jet');
    opt.viewAngle = getOpt(options, 'viewAngle', [135, 30]);
    opt.lighting  = getOpt(options, 'lighting', true);
    opt.showSTL   = getOpt(options, 'showSTL', false);
    opt.stlModel  = getOpt(options, 'stlModel', '');
    opt.saveFigs  = getOpt(options, 'saveFigs', false);
    opt.saveDir   = getOpt(options, 'saveDir', fullfile('..', 'results'));
    opt.figPrefix = getOpt(options, 'figPrefix', data.fileName);
end

function val = getOpt(s, field, defaultVal)
    if isfield(s, field) && ~isempty(s.(field))
        val = s.(field);
    else
        val = defaultVal;
    end
end

%% ========================================================================
%  PLOTSPHERICAL  3D spherical polar RCS surface
% ========================================================================
function plotSpherical(data, rcsData, compLabel, cIdx, nComp, opt)
    figName = sprintf('FFE - 3D Spherical [%s] %s', data.fileName, compLabel);
    fig = figure('Name', figName, 'NumberTitle', 'off', ...
                 'Position', [80 + (cIdx-1)*30, 60, 750, 650], 'Color', 'w');

    % Get angular coordinates
    thVec = deg2rad(data.theta(1, :));
    phVec = deg2rad(data.phi(:, 1));

    % Close the phi wrap for a seamless sphere
    phFull  = [phVec; phVec(1) + 2*pi];
    rcsFull = rcsData;
    rcsFull(end+1, :) = rcsData(1, :);

    [ThetaGrid, PhiGrid] = meshgrid(thVec, phFull);

    % Scale RCS to positive range for spherical radius mapping
    minRCS = min(rcsFull(:));
    rMin   = 0.3;
    rScale = rcsFull - minRCS + rMin;

    % Convert to Cartesian for surf
    X = rScale .* sin(ThetaGrid) .* cos(PhiGrid);
    Y = rScale .* sin(ThetaGrid) .* sin(PhiGrid);
    Z = rScale .* cos(ThetaGrid);

    % Surface plot
    surf(X, Y, Z, rcsFull, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
    hold on;

    % Reference sphere at rMin
    [Sx, Sy, Sz] = sphere(40);
    Sx = Sx * rMin; Sy = Sy * rMin; Sz = Sz * rMin;
    mesh(Sx, Sy, Sz, 'FaceColor', 'none', 'EdgeColor', [0.6 0.6 0.6], ...
         'EdgeAlpha', 0.2, 'LineWidth', 0.4);

    % Coordinate axes
    mx = max(abs([X(:); Y(:); Z(:)])) * 1.2;
    line([-mx mx], [0 0], [0 0], 'Color', 'k', 'LineWidth', 0.5);
    line([0 0], [-mx mx], [0 0], 'Color', 'k', 'LineWidth', 0.5);
    line([0 0], [0 0], [-mx mx], 'Color', 'k', 'LineWidth', 0.5);
    text(mx*1.08, 0, 0, 'X', 'FontSize', 11, 'FontWeight', 'bold');
    text(0, mx*1.08, 0, 'Y', 'FontSize', 11, 'FontWeight', 'bold');
    text(0, 0, mx*1.08, 'Z', 'FontSize', 11, 'FontWeight', 'bold');

    % Styling
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal; axis tight;
    view(opt.viewAngle);
    if opt.lighting
        lighting gouraud;
        camlight('headlight');
        material dull;
    end
    colormap(opt.colormap);
    cb = colorbar;
    cb.Label.String = 'RCS (dBsm)';
    cb.Label.FontSize = 11;

    % Title with metadata
    titleStr = sprintf('3D Spherical RCS  %s  [%s]', compLabel, data.fileName);
    if isfield(data, 'paramStruct')
        ps = data.paramStruct;
        if isfield(ps, 'Mode'), titleStr = [titleStr '  |  ' ps.Mode]; end %#ok<AGROW>
    end
    title(titleStr, 'FontSize', 12, 'Interpreter', 'tex');

    if opt.saveFigs
        saveFigure(fig, fullfile(opt.saveDir, ...
                   sprintf('%s_rcs_3d_sphere_%s', opt.figPrefix, compLabel)), opt);
    end
end

%% ========================================================================
%  PLOTTHETAPHIMESH  Theta-Phi 3D mesh/surface (waterfall style)
% ========================================================================
function plotThetaPhiMesh(data, rcsData, compLabel, cIdx, nComp, opt)
    figName = sprintf('FFE - 3D Mesh [%s] %s', data.fileName, compLabel);
    fig = figure('Name', figName, 'NumberTitle', 'off', ...
                 'Position', [80 + (cIdx-1)*30, 40, 820, 580], 'Color', 'w');

    % Prepare meshgrid
    phiDeg   = data.phi;
    thetaDeg = data.theta;
    plotData = rcsData;

    % Wrap-around for full 360° display if partial scan
    phiRange   = max(phiDeg(:)) - min(phiDeg(:));
    thetaRange = max(thetaDeg(:)) - min(thetaDeg(:));
    if phiRange < 350 && data.ip > 1
        phiDeg   = [phiDeg; phiDeg(1, :) + 360];
        plotData = [plotData; plotData(1, :)];
    end

    [T, P] = meshgrid(thetaDeg(1, :), phiDeg(:, 1));
    if size(plotData, 1) ~= size(P, 1)
        plotData = plotData';
    end

    % Surface plot
    surf(T, P, plotData, 'EdgeColor', 'interp', 'FaceColor', 'interp');
    hold on;

    % Peak marker (draw a vertical line from floor to peak)
    [maxVal, maxIdx] = max(plotData(:));
    [maxRow, maxCol] = ind2sub(size(plotData), maxIdx);
    plot3(T(maxRow, maxCol), P(maxRow, maxCol), maxVal, ...
          'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineWidth', 1.5);

    % Styling
    xlabel('\theta (deg)', 'FontSize', 12);
    ylabel('\phi (deg)', 'FontSize', 12);
    zlabel('RCS (dBsm)', 'FontSize', 12);
    title(sprintf('RCS Theta-Phi Pattern  %s  [%s]', compLabel, data.fileName), ...
          'FontSize', 12, 'Interpreter', 'tex');
    view(opt.viewAngle);
    colormap(opt.colormap);
    cb = colorbar;
    cb.Label.String = 'RCS (dBsm)';
    cb.Label.FontSize = 11;
    if opt.lighting
        lighting gouraud;
        camlight('headlight');
    end
    grid on;

    if opt.saveFigs
        saveFigure(fig, fullfile(opt.saveDir, ...
                   sprintf('%s_rcs_3d_mesh_%s', opt.figPrefix, compLabel)), opt);
    end
end

%% ========================================================================
%  PLOTHEATMAP  2D pcolor heatmap with contour overlay
% ========================================================================
function plotHeatmap(data, rcsData, compLabel, cIdx, nComp, opt)
    figName = sprintf('FFE - Heatmap [%s] %s', data.fileName, compLabel);
    fig = figure('Name', figName, 'NumberTitle', 'off', ...
                 'Position', [80 + (cIdx-1)*30, 50, 820, 550], 'Color', 'w');

    if data.is2D
        % ---- 2D heatmap ----
        subplot(1, 2, 1);
        thetaDeg = data.theta(1, :);
        phiDeg   = data.phi(:, 1);
        [T, P] = meshgrid(thetaDeg, phiDeg);

        pcolor(T, P, rcsData);
        shading interp;
        hold on;

        % Contour overlay
        nContours = 10;
        contour(T, P, rcsData, nContours, 'k-', 'LineWidth', 0.5);

        % Peak marker
        [maxVal, maxIdx] = max(rcsData(:));
        [maxRow, maxCol] = ind2sub(size(rcsData), maxIdx);
        plot(T(maxRow, maxCol), P(maxRow, maxCol), 'wx', ...
             'MarkerSize', 10, 'LineWidth', 2);
        text(T(maxRow, maxCol), P(maxRow, maxCol), ...
             sprintf('  %.1f dBsm', maxVal), 'Color', 'w', ...
             'FontSize', 10, 'FontWeight', 'bold');

        xlabel('\theta (deg)', 'FontSize', 12);
        ylabel('\phi (deg)', 'FontSize', 12);
        title(sprintf('Heatmap + Contour  %s', compLabel), 'FontSize', 11);
        colormap(gca, opt.colormap);
        cb1 = colorbar;
        cb1.Label.String = 'RCS (dBsm)';
        axis tight; grid on;

        % ---- 3D Waterfall (tight) ----
        subplot(1, 2, 2);
        waterfall(T, P, rcsData);
        xlabel('\theta (deg)', 'FontSize', 12);
        ylabel('\phi (deg)', 'FontSize', 12);
        zlabel('RCS (dBsm)', 'FontSize', 12);
        title(sprintf('Waterfall  %s', compLabel), 'FontSize', 11);
        view(opt.viewAngle);
        colormap(gca, opt.colormap);
        grid on;

        sgtitle(sprintf('RCS Pattern Analysis  —  %s  |  %s', data.fileName, compLabel), ...
                'FontSize', 13);

    else
        % ---- 1D profile as filled area plot ----
        if data.ip == 1
            angle = data.theta(1, :);
            xLab = '\theta (deg)';
        else
            angle = data.phi(:, 1)';
            xLab = '\phi (deg)';
        end

        plot(angle, rcsData, '-', 'LineWidth', 2.5, 'Color', [0.2 0.4 0.8]);
        hold on;
        fill([angle(1), angle, angle(end)], ...
             [min(rcsData)-5, rcsData, min(rcsData)-5], ...
             [0.2 0.4 0.8], 'FaceAlpha', 0.15, 'EdgeColor', 'none');

        % Peak marker
        [maxVal, maxIdx] = max(rcsData);
        plot(angle(maxIdx), maxVal, 'ro', ...
             'MarkerSize', 10, 'MarkerFaceColor', 'r');
        text(angle(maxIdx), maxVal, ...
             sprintf('  %.1f dBsm @ %.0f°', maxVal, angle(maxIdx)), ...
             'FontSize', 10, 'FontWeight', 'bold');

        xlabel(xLab, 'FontSize', 13);
        ylabel('RCS (dBsm)', 'FontSize', 13);
        title(sprintf('RCS Profile  %s  [%s]', compLabel, data.fileName), 'FontSize', 12);
        grid on; box on;
        xlim([min(angle) max(angle)]);
    end

    if opt.saveFigs
        saveFigure(fig, fullfile(opt.saveDir, ...
                   sprintf('%s_rcs_heatmap_%s', opt.figPrefix, compLabel)), opt);
    end
end

%% ========================================================================
%  SAVEFIGURE  Save figure to file
% ========================================================================
function saveFigure(fig, basePath, opt)
    if ~exist(opt.saveDir, 'dir')
        [status, msg] = mkdir(opt.saveDir);
        if ~status
            warning('ffe_plot_rcs_3d:CannotCreateDir', ...
                    'Cannot create save directory: %s', msg);
            return;
        end
    end

    try
        exportgraphics(fig, [basePath '.png'], 'Resolution', 300);
        fprintf('  Figure saved: %s.png\n', basePath);
    catch
        try
            saveas(fig, [basePath '.png']);
            fprintf('  Figure saved: %s.png\n', basePath);
        catch ME
            warning('ffe_plot_rcs_3d:SaveFailed', ...
                    'Failed to save figure: %s', ME.message);
        end
    end
end
