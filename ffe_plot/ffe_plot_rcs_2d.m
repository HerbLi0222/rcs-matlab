% FFE_PLOT_RCS_2D  Enhanced 2D RCS pattern visualization
%
%   Creates publication-quality 2D RCS plots for 1D angular cuts:
%     - Cartesian plot with peak markers, -3dB beamwidth annotation
%     - Polar plot with dB-scale rings and main lobe indicator
%
%   For 2D grid data, extracts and plots the center-cut slices.
%
%   Input:
%       data    - data struct from ffe_load_data
%       options - (optional) struct with fields:
%           .component   - 'theta' (default), 'phi', or 'both'
%           .smoothWin   - smoothing window (0 = no smoothing, default = 0)
%           .polarPlot   - true (default) to generate polar plot
%           .cartesianPlot - true (default) to generate cartesian plot
%           .showBeamwidth - true (default) to annotate -3dB beamwidth
%           .showPeaks   - true (default) to mark local peaks
%           .peakThreshold - prominence threshold for peak detection (dB, default = 3)
%           .saveFigs    - true/false to save figures to file (default = false)
%           .saveDir     - directory for saved figures (default = 'results/')
%           .figPrefix   - prefix for saved figure filenames
%
%   Usage:
%     >> data = ffe_load_data('results/temp_20260616220000.dat');
%     >> ffe_plot_rcs_2d(data);
%     >> ffe_plot_rcs_2d(data, struct('component', 'both', 'smoothWin', 3));
%
%   See also: ffe_main, ffe_load_data, ffe_plot_rcs_3d, ffe_utils

function ffe_plot_rcs_2d(data, options)

    %% ---- Default options ----
    if nargin < 2, options = struct(); end
    opt = setDefaults(options, data);

    %% ---- Validate data ----
    if isempty(data.Sth)
        error('ffe_plot_rcs_2d:NoData', ...
              'No RCS data available for plotting. Load a .dat or .mat file first.');
    end

    %% ---- Extract 1D slices if 2D grid ----
    if data.is2D
        fprintf('  Data is 2D grid — extracting center cuts for 1D plotting.\n');
        [angle, rcsTh, rcsPh, angLabel, sliceLabel] = extract1DCut(data);
    else
        [angle, rcsTh, rcsPh, angLabel, sliceLabel] = extract1DCut(data);
    end

    %% ---- Apply smoothing if requested ----
    if opt.smoothWin > 0
        rcsTh = smooth_rcs(rcsTh, opt.smoothWin);
        rcsPh = smooth_rcs(rcsPh, opt.smoothWin);
    end

    %% ---- Cartesian plot ----
    if opt.cartesianPlot
        plotCartesian(angle, rcsTh, rcsPh, angLabel, sliceLabel, data, opt);
    end

    %% ---- Polar plot ----
    if opt.polarPlot
        plotPolar(angle, rcsTh, rcsPh, angLabel, sliceLabel, data, opt);
    end
end

%% ========================================================================
%  SETDEFAULTS  Fill in default option values
% ========================================================================
function opt = setDefaults(options, data)
    opt.component      = getOpt(options, 'component', 'theta');
    opt.smoothWin      = getOpt(options, 'smoothWin', 0);
    opt.polarPlot      = getOpt(options, 'polarPlot', true);
    opt.cartesianPlot  = getOpt(options, 'cartesianPlot', true);
    opt.showBeamwidth  = getOpt(options, 'showBeamwidth', true);
    opt.showPeaks      = getOpt(options, 'showPeaks', true);
    opt.peakThreshold  = getOpt(options, 'peakThreshold', 3);
    opt.saveFigs       = getOpt(options, 'saveFigs', false);
    opt.saveDir        = getOpt(options, 'saveDir', fullfile('..', 'results'));
    opt.figPrefix      = getOpt(options, 'figPrefix', data.fileName);
end

function val = getOpt(s, field, defaultVal)
    if isfield(s, field) && ~isempty(s.(field))
        val = s.(field);
    else
        val = defaultVal;
    end
end

%% ========================================================================
%  EXTRACT1DCUT  Extract 1D angular cut from data
% ========================================================================
function [angle, rcsTh, rcsPh, angLabel, sliceLabel] = extract1DCut(data)
    if data.ip == 1
        % Single phi, sweep theta
        angle     = data.theta(1, :);
        rcsTh     = data.Sth(1, :);
        if size(data.Sph, 1) >= 1
            rcsPh = data.Sph(1, :);
        else
            rcsPh = zeros(size(rcsTh));
        end
        angLabel   = '\theta (deg)';
        sliceLabel = sprintf('\\phi = %.0f°', data.phi(1, 1));
    elseif data.it == 1
        % Single theta, sweep phi
        angle     = data.phi(:, 1)';
        rcsTh     = data.Sth(:, 1)';
        rcsPh     = data.Sph(:, 1)';
        angLabel   = '\phi (deg)';
        sliceLabel = sprintf('\\theta = %.0f°', data.theta(1, 1));
    else
        % 2D grid: take center theta cut at center phi
        midPhi = round(data.ip / 2);
        angle  = data.theta(midPhi, :);
        rcsTh  = data.Sth(midPhi, :);
        rcsPh  = data.Sph(midPhi, :);
        angLabel   = '\theta (deg)';
        sliceLabel = sprintf('\\phi = %.0f° (center cut)', data.phi(midPhi, 1));
    end

    % Ensure row vectors
    angle = angle(:)';
    rcsTh = rcsTh(:)';
    rcsPh = rcsPh(:)';
end

%% ========================================================================
%  PLOTCARTESIAN  Generate enhanced cartesian RCS plot
% ========================================================================
function plotCartesian(angle, rcsTh, rcsPh, angLabel, sliceLabel, data, opt)
    fig = setup_figure_custom( ...
        sprintf('FFE - RCS Cartesian  [%s]', data.fileName), 850, 500);

    plotComponent = lower(opt.component);

    hold on;
    legends = {};
    colors = lines(4);

    % ---- Plot S_theta ----
    if ismember(plotComponent, {'theta', 'both'})
        h1 = plot(angle, rcsTh, '-', 'LineWidth', 2.0, 'Color', colors(1,:));
        legends{end+1} = 'S_{\theta}';

        % Annotate peaks on S_theta
        if opt.showPeaks && length(angle) > 3
            [peakVals, peakIdx, peakProm] = find_peaks_rcs(rcsTh, opt.peakThreshold);
            for p = 1:length(peakIdx)
                plot(angle(peakIdx(p)), peakVals(p), 'v', ...
                     'MarkerSize', 8, 'MarkerFaceColor', colors(1,:), ...
                     'MarkerEdgeColor', 'k');
            end
            if ~isempty(peakIdx)
                [~, mainIdx] = max(peakVals);
                text(angle(peakIdx(mainIdx)), peakVals(mainIdx), ...
                     sprintf('  %.1f dBsm @ %.0f°', peakVals(mainIdx), angle(peakIdx(mainIdx))), ...
                     'FontSize', 9, 'Color', colors(1,:), 'FontWeight', 'bold');
            end
        end

        % Annotate -3dB beamwidth
        if opt.showBeamwidth && length(angle) > 3
            [bw, leftI, rightI] = beamwidth(rcsTh, angle, 3);
            if ~isnan(bw)
                [peakVal, ~] = max(rcsTh);
                level3dB = peakVal - 3;
                yLim = ylim;
                % Draw horizontal line at -3dB level
                plot([angle(1) angle(end)], [level3dB level3dB], '--', ...
                     'Color', [0.5 0.5 0.5], 'LineWidth', 1);
                % Draw vertical markers at beam edges
                if ~isnan(leftI)
                    leftAng = interp1(1:length(angle), angle, leftI, 'linear', 'extrap');
                    plot([leftAng leftAng], yLim, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
                end
                if ~isnan(rightI)
                    rightAng = interp1(1:length(angle), angle, rightI, 'linear', 'extrap');
                    plot([rightAng rightAng], yLim, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
                end
                % Beamwidth annotation
                midAng = (angle(1) + angle(end)) / 2;
                text(midAng, level3dB - 1.5, ...
                     sprintf('-3dB BW: %.1f°', bw), ...
                     'FontSize', 9, 'Color', [0.4 0.4 0.4], ...
                     'HorizontalAlignment', 'center');
            end
        end
    end

    % ---- Plot S_phi ----
    if ismember(plotComponent, {'phi', 'both'})
        h2 = plot(angle, rcsPh, '--', 'LineWidth', 2.0, 'Color', colors(2,:));
        legends{end+1} = 'S_{\phi}';

        if opt.showPeaks && length(angle) > 3
            [peakVals, peakIdx] = find_peaks_rcs(rcsPh, opt.peakThreshold);
            for p = 1:length(peakIdx)
                plot(angle(peakIdx(p)), peakVals(p), '^', ...
                     'MarkerSize', 6, 'MarkerFaceColor', colors(2,:), ...
                     'MarkerEdgeColor', 'k');
            end
        end
    end

    % ---- Labels and styling ----
    xlabel(angLabel, 'FontSize', 13);
    ylabel('RCS (dBsm)', 'FontSize', 13);

    % Build title with metadata
    titleStr = sprintf('RCS vs %s  [%s]', angLabel, sliceLabel);
    if isfield(data, 'paramStruct') && ~isempty(fieldnames(data.paramStruct))
        ps = data.paramStruct;
        % Try to append mode/freq/polarization if available
        if isfield(ps, 'Mode'), titleStr = [titleStr sprintf('  |  %s', ps.Mode)]; end %#ok<AGROW>
        if isfield(ps, 'RadarFrequency')
            titleStr = [titleStr sprintf('  |  %s GHz', ps.RadarFrequency)]; %#ok<AGROW>
        end
    end
    title(titleStr, 'FontSize', 12, 'Interpreter', 'tex');

    legend(legends, 'Location', 'best', 'FontSize', 11);

    xlim([min(angle) max(angle)]);
    grid on; box on;

    % ---- Save if requested ----
    if opt.saveFigs
        saveFigure(fig, fullfile(opt.saveDir, ...
                   [opt.figPrefix '_rcs_2d_cartesian']), opt);
    end
end

%% ========================================================================
%  PLOTPOLAR  Generate enhanced polar RCS plot
% ========================================================================
function plotPolar(angle, rcsTh, rcsPh, angLabel, sliceLabel, data, opt)
    % Create figure and let polarplot create its own PolarAxes
    fig = figure('Name', sprintf('FFE - RCS Polar  [%s]', data.fileName), ...
                 'NumberTitle', 'off', ...
                 'Position', [100, 100, 650, 580], ...
                 'Color', 'w');

    plotComponent = lower(opt.component);

    % Normalize data for polar display (shift to positive range)
    allVals = [];
    if ismember(plotComponent, {'theta', 'both'}), allVals = [allVals, rcsTh]; end
    if ismember(plotComponent, {'phi', 'both'}),   allVals = [allVals, rcsPh]; end
    minAll = min(allVals);
    offset = abs(minAll) + 5;  % shift so minimum is above 0 with margin

    legends = {};
    colors = lines(4);

    % ---- First polarplot call creates the PolarAxes ----
    if ismember(plotComponent, {'theta', 'both'})
        rTh = rcsTh + offset;
        polarplot(deg2rad(angle), rTh, '-', 'LineWidth', 2.0, 'Color', colors(1,:));
        legends{end+1} = 'S_{\theta}';
        hold on;

        % Mark main lobe peak
        [peakVal, peakIdx] = max(rcsTh);
        polarplot(deg2rad(angle(peakIdx)), rTh(peakIdx), 'v', ...
                  'MarkerSize', 8, 'MarkerFaceColor', colors(1,:), ...
                  'MarkerEdgeColor', 'k');
    end

    % ---- Plot S_phi on polar axes ----
    if ismember(plotComponent, {'phi', 'both'})
        rPh = rcsPh + offset;
        polarplot(deg2rad(angle), rPh, '--', 'LineWidth', 2.0, 'Color', colors(2,:));
        legends{end+1} = 'S_{\phi}';
        if isempty(legends) || length(legends) <= 1
            hold on;
        end
    end

    % ---- Title and legend ----
    titleStr = sprintf('Polar RCS Pattern  [%s]  |  dB offset = %.0f', sliceLabel, offset);
    title(titleStr, 'FontSize', 12, 'Interpreter', 'tex');
    if ~isempty(legends)
        legend(legends, 'Location', 'bestoutside', 'FontSize', 10);
    end

    % ---- Save if requested ----
    if opt.saveFigs
        saveFigure(fig, fullfile(opt.saveDir, ...
                   [opt.figPrefix '_rcs_2d_polar']), opt);
    end
end

%% ========================================================================
%  SETUP_FIGURE_CUSTOM  Create a figure window (local version for standalone use)
% ========================================================================
function fig = setup_figure_custom(name, width, height)
    fig = figure('Name', name, 'NumberTitle', 'off', ...
                 'Position', [100, 100, width, height], ...
                 'Color', 'w');
    set(gca, 'FontSize', 11, 'LineWidth', 1, 'Box', 'on');
    grid on;
    hold on;
end

%% ========================================================================
%  SAVEFIGURE  Save figure to file
% ========================================================================
function saveFigure(fig, basePath, opt)
    if ~exist(opt.saveDir, 'dir')
        [status, msg] = mkdir(opt.saveDir);
        if ~status
            warning('ffe_plot_rcs_2d:CannotCreateDir', ...
                    'Cannot create save directory: %s', msg);
            return;
        end
    end

    % Try exportgraphics (R2020a+), fall back to saveas
    try
        exportgraphics(fig, [basePath '.png'], 'Resolution', 300);
        fprintf('  Figure saved: %s.png\n', basePath);
    catch
        try
            saveas(fig, [basePath '.png']);
            fprintf('  Figure saved: %s.png\n', basePath);
        catch ME
            warning('ffe_plot_rcs_2d:SaveFailed', ...
                    'Failed to save figure: %s', ME.message);
        end
    end
end
