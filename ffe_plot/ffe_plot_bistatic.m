% FFE_PLOT_BISTATIC  Bistatic RCS visualization for POSTFEKO .ffe data
%   双站RCS方向图可视化
%
%   Creates comprehensive bistatic RCS visualizations from multi-block
%   POSTFEKO .ffe data:
%     - Bistatic RCS heatmap (observation angle vs incident angle)
%     - Monostatic cut (theta_obs = theta_inc)
%     - Fixed-incident-angle far-field cuts
%     - Forward/backward scattering comparison
%     - 3D surface of bistatic RCS matrix
%
%   Input:
%       data    - data struct from ffe_load_data (must be .ffe / bistatic)
%       options - (optional) struct with fields:
%           .incidentAngles - incident angle indices to plot [array] (default = [1, center, end])
%           .plotType       - 'all' (default), 'heatmap', 'monostatic', 'cuts', 'forward_back', '3d'
%           .component      - 'total' (default), 'theta', 'phi'
%           .colormap       - 'jet' (default), 'parula', 'hot', 'turbo'
%           .saveFigs       - true/false to save figures (default = false)
%           .saveDir        - directory for saved figures
%           .figPrefix      - prefix for saved figure filenames
%
%   Usage:
%     >> data = ffe_load_data('ffe_file/result_rocket_po.ffe');
%     >> ffe_plot_bistatic(data);
%     >> ffe_plot_bistatic(data, struct('plotType', 'monostatic'));
%
%   See also: ffe_main, ffe_load_data, ffe_plot_rcs_2d

function ffe_plot_bistatic(data, options)

    %% ---- Validate ----
    if ~isfield(data, 'bistatic') || ~data.bistatic
        error('ffe_plot_bistatic:NotBistatic', ...
              'This function requires bistatic .ffe data. Load a .ffe file first.');
    end

    if nargin < 2, options = struct(); end
    opt = setDefaults(options, data);

    %% ---- Extract data for selected component ----
    switch lower(opt.component)
        case 'theta'
            rcsLin = data.RCS_theta;
            compLabel = 'S_{\theta}';
        case 'phi'
            rcsLin = data.RCS_phi;
            compLabel = 'S_{\phi}';
        case 'total'
            rcsLin = data.RCS_total;
            compLabel = 'S_{total}';
    end

    rcsDB = 10 * log10(max(rcsLin, 1e-30));

    thetaObs = data.theta(1, :)';  % observation theta [deg] (it x 1)
    thetaInc = data.theta_inc(:);  % incident theta [deg] (N_inc x 1)
    it = data.it;
    N_inc = data.N_inc;

    %% ---- Figure 1: Bistatic Heatmap ----
    if ismember(opt.plotType, {'all', 'heatmap'})
        plotBistaticHeatmap(thetaObs, thetaInc, rcsDB, it, N_inc, ...
                            compLabel, data, opt);
    end

    %% ---- Figure 2: Monostatic Cut ----
    if ismember(opt.plotType, {'all', 'monostatic'})
        plotMonostaticCut(thetaObs, thetaInc, rcsDB, compLabel, data, opt);
    end

    %% ---- Figure 3: Fixed-Incident-Angle Cuts ----
    if ismember(opt.plotType, {'all', 'cuts'})
        plotIncidentCuts(thetaObs, thetaInc, rcsDB, compLabel, data, opt);
    end

    %% ---- Figure 4: Forward vs Backward Scattering ----
    if ismember(opt.plotType, {'all', 'forward_back'})
        plotForwardBackward(thetaObs, thetaInc, rcsDB, compLabel, data, opt);
    end

    %% ---- Figure 5: 3D Surface ----
    if ismember(opt.plotType, {'all', '3d'})
        plotBistatic3D(thetaObs, thetaInc, rcsDB, compLabel, data, opt);
    end
end

%% ========================================================================
function opt = setDefaults(options, data)
    opt.component = getOpt(options, 'component', 'total');
    opt.plotType  = getOpt(options, 'plotType', 'all');
    opt.colormap  = getOpt(options, 'colormap', 'jet');
    opt.saveFigs  = getOpt(options, 'saveFigs', false);
    opt.saveDir   = getOpt(options, 'saveDir', './');
    opt.figPrefix = getOpt(options, 'figPrefix', data.fileName);

    % Default incident angle indices for cuts plot
    N = data.N_inc;
    opt.incidentAngles = getOpt(options, 'incidentAngles', ...
        unique([1, round(N/4), round(N/2), round(3*N/4), N]));
    % Ensure valid indices
    opt.incidentAngles = opt.incidentAngles(opt.incidentAngles >= 1 & opt.incidentAngles <= N);
end

function val = getOpt(s, field, defaultVal)
    if isfield(s, field) && ~isempty(s.(field))
        val = s.(field);
    else
        val = defaultVal;
    end
end

%% ========================================================================
%  FIGURE 1: Bistatic RCS Heatmap (observation angle vs incident angle)
% ========================================================================
function plotBistaticHeatmap(thetaObs, thetaInc, rcsDB, it, N_inc, ...
                             compLabel, data, opt)
    fig = figure('Name', sprintf('FFE - Bistatic Heatmap [%s]', data.fileName), ...
                 'NumberTitle', 'off', ...
                 'Position', [50, 80, 950, 600], 'Color', 'w');

    % ---- Main heatmap ----
    subplot(1, 2, 1);
    imagesc(thetaInc, thetaObs, rcsDB);
    set(gca, 'YDir', 'normal');
    xlabel('Incident Angle \theta_i (deg)', 'FontSize', 12);
    ylabel('Observation Angle \theta_o (deg)', 'FontSize', 12);
    title(sprintf('Bistatic RCS  %s  —  %s', compLabel, data.fileName), ...
          'FontSize', 11, 'Interpreter', 'none');
    colormap(gca, opt.colormap);
    cb1 = colorbar;
    cb1.Label.String = 'RCS (dBsm)';
    cb1.Label.FontSize = 11;
    axis tight;

    % Monostatic line (theta_obs = theta_inc)
    hold on;
    monoMin = max(min(thetaInc), min(thetaObs));
    monoMax = min(max(thetaInc), max(thetaObs));
    plot([monoMin monoMax], [monoMin monoMax], 'k--', 'LineWidth', 1.5);
    text(monoMax, monoMax, ' Mono', 'FontSize', 9, 'Color', 'k', ...
         'VerticalAlignment', 'bottom');

    % Forward scatter line (theta_obs = theta_inc, phi_obs = 0)
    % For 2D: shows specular direction

    % ---- Diagonal profile (monostatic cut) ----
    subplot(1, 2, 2);
    [monoAngles, monoRCS] = extractMonostatic(thetaObs, thetaInc, rcsDB);
    plot(monoAngles, monoRCS, 'r-', 'LineWidth', 2.5);
    xlabel('Angle \theta (deg)', 'FontSize', 12);
    ylabel('Monostatic RCS (dBsm)', 'FontSize', 12);
    title(sprintf('Monostatic Cut (\\theta_o = \\theta_i)  %s', compLabel), ...
          'FontSize', 11);
    grid on; box on;

    % Peak marker
    [peakVal, peakIdx] = max(monoRCS);
    hold on;
    plot(monoAngles(peakIdx), peakVal, 'bo', ...
         'MarkerSize', 8, 'MarkerFaceColor', 'b');
    text(monoAngles(peakIdx), peakVal, ...
         sprintf('  %.1f dBsm @ %.0f°', peakVal, monoAngles(peakIdx)), ...
         'FontSize', 9, 'Color', 'b');

    sgtitle(sprintf('Bistatic RCS Analysis  —  %s  |  f = %.2f GHz', ...
            data.fileName, data.freq_ghz), 'FontSize', 12, 'Interpreter', 'none');

    if opt.saveFigs
        saveFig(fig, fullfile(opt.saveDir, [opt.figPrefix '_bistatic_heatmap']), opt);
    end
end

%% ========================================================================
%  FIGURE 2: Monostatic Cut (diagonal of bistatic matrix)
% ========================================================================
function plotMonostaticCut(thetaObs, thetaInc, rcsDB, compLabel, data, opt)
    fig = figure('Name', sprintf('FFE - Monostatic Cut [%s]', data.fileName), ...
                 'NumberTitle', 'off', ...
                 'Position', [80, 100, 800, 500], 'Color', 'w');

    [monoAngles, monoRCS] = extractMonostatic(thetaObs, thetaInc, rcsDB);

    % ---- Monostatic RCS plot ----
    subplot(1, 2, 1);
    plot(monoAngles, monoRCS, 'r-', 'LineWidth', 2.5);
    xlabel('\theta (deg)', 'FontSize', 12);
    ylabel('Monostatic RCS (dBsm)', 'FontSize', 12);
    title(sprintf('Monostatic RCS  %s', compLabel), 'FontSize', 12);
    grid on; box on;
    xlim([min(monoAngles) max(monoAngles)]);

    % Peak annotation
    [peakVal, peakIdx] = max(monoRCS);
    hold on;
    plot(monoAngles(peakIdx), peakVal, 'bo', ...
         'MarkerSize', 10, 'MarkerFaceColor', 'b');
    text(monoAngles(peakIdx), peakVal, ...
         sprintf('  Peak: %.1f dBsm @ %.0f°', peakVal, monoAngles(peakIdx)), ...
         'FontSize', 10, 'FontWeight', 'bold', 'Color', 'b');

    % ---- Polar plot ----
    subplot(1, 2, 2);
    minRCS = min(monoRCS);
    rPolar = monoRCS - minRCS + 5;
    polarplot(deg2rad(monoAngles), rPolar, 'r-', 'LineWidth', 2);
    hold on;
    [peakVal, peakIdx] = max(monoRCS);
    polarplot(deg2rad(monoAngles(peakIdx)), rPolar(peakIdx), 'bo', ...
              'MarkerSize', 8, 'MarkerFaceColor', 'b');
    title(sprintf('Monostatic Polar  %s', compLabel), 'FontSize', 12);

    sgtitle(sprintf('Monostatic RCS (\\theta_o = \\theta_i)  —  %s  |  f = %.2f GHz', ...
            data.fileName, data.freq_ghz), 'FontSize', 12, 'Interpreter', 'none');

    if opt.saveFigs
        saveFig(fig, fullfile(opt.saveDir, [opt.figPrefix '_monostatic']), opt);
    end
end

%% ========================================================================
%  FIGURE 3: Fixed-Incident-Angle Far-Field Cuts
% ========================================================================
function plotIncidentCuts(thetaObs, thetaInc, rcsDB, compLabel, data, opt)
    fig = figure('Name', sprintf('FFE - Incident Cuts [%s]', data.fileName), ...
                 'NumberTitle', 'off', ...
                 'Position', [40, 60, 1100, 550], 'Color', 'w');

    idxList = opt.incidentAngles;
    N_cuts = length(idxList);
    colors = lines(N_cuts);

    % ---- Overlay plot ----
    subplot(1, 2, 1);
    hold on;
    legends = {};
    for i = 1:N_cuts
        idx = idxList(i);
        rcsCut = rcsDB(:, idx);
        if size(rcsCut, 1) ~= length(thetaObs)
            rcsCut = rcsCut';
        end
        plot(thetaObs, rcsCut, '-', 'LineWidth', 1.8, 'Color', colors(i,:));
        legends{end+1} = sprintf('\\theta_i = %.0f°', thetaInc(idx)); %#ok<AGROW>
    end
    xlabel('Observation Angle \theta_o (deg)', 'FontSize', 12);
    ylabel('RCS (dBsm)', 'FontSize', 12);
    title(sprintf('Far-Field Cuts at Fixed Incident Angles  %s', compLabel), ...
          'FontSize', 11);
    legend(legends, 'Location', 'bestoutside', 'FontSize', 9);
    grid on; box on;
    xlim([min(thetaObs) max(thetaObs)]);

    % ---- Waterfall (offset for visibility) ----
    subplot(1, 2, 2);
    hold on;
    offset = 0;
    offsetStep = (max(rcsDB(:)) - min(rcsDB(:))) * 0.15;
    for i = 1:N_cuts
        idx = idxList(i);
        rcsCut = rcsDB(:, idx);
        if size(rcsCut, 1) ~= length(thetaObs)
            rcsCut = rcsCut';
        end
        plot3(thetaObs, repmat(thetaInc(idx), size(thetaObs)), rcsCut, ...
              '-', 'LineWidth', 1.8, 'Color', colors(i,:));
    end
    xlabel('\theta_o (deg)', 'FontSize', 12);
    ylabel('\theta_i (deg)', 'FontSize', 12);
    zlabel('RCS (dBsm)', 'FontSize', 12);
    title('Waterfall View', 'FontSize', 11);
    view([-45, 30]);
    grid on;

    sgtitle(sprintf('Far-Field Cuts vs Incident Angle  —  %s  |  f = %.2f GHz', ...
            data.fileName, data.freq_ghz), 'FontSize', 12, 'Interpreter', 'none');

    if opt.saveFigs
        saveFig(fig, fullfile(opt.saveDir, [opt.figPrefix '_incident_cuts']), opt);
    end
end

%% ========================================================================
%  FIGURE 4: Forward vs Backward Scattering
% ========================================================================
function plotForwardBackward(thetaObs, thetaInc, rcsDB, compLabel, data, opt)
    fig = figure('Name', sprintf('FFE - Forward/Backward [%s]', data.fileName), ...
                 'NumberTitle', 'off', ...
                 'Position', [50, 80, 1000, 520], 'Color', 'w');

    % Forward scatter: observation near incident direction
    % Backward scatter: observation opposite to incident direction

    % For each incident angle, find RCS at specular (forward) and monostatic (backward)
    N_inc = data.N_inc;
    forwardRCS  = zeros(N_inc, 1);
    backwardRCS = zeros(N_inc, 1);
    forwardAng  = zeros(N_inc, 1);

    for i = 1:N_inc
        th_i = thetaInc(i);

        % Monostatic (backward): theta_obs = theta_inc
        [~, monoIdx] = min(abs(thetaObs - th_i));
        if monoIdx <= length(thetaObs) && monoIdx <= size(rcsDB, 1)
            backwardRCS(i) = rcsDB(monoIdx, i);
        end

        % Forward scatter: theta_obs ≈ theta_inc (same direction, transmission through target)
        % In far-field, forward scatter is at theta_obs = theta_inc
        forwardRCS(i) = backwardRCS(i);  % same as mono in 1D

        % Bistatic forward: theta_obs = -theta_inc (opposite direction)
        [~, fwdIdx] = min(abs(thetaObs + th_i));
        if fwdIdx <= length(thetaObs) && fwdIdx <= size(rcsDB, 1)
            forwardAng(i) = thetaObs(fwdIdx);
        end
    end

    % ---- Monostatic (backscatter) vs incident angle ----
    subplot(1, 2, 1);
    plot(thetaInc, backwardRCS, 'b-', 'LineWidth', 2.5);
    xlabel('Incident Angle \theta_i (deg)', 'FontSize', 12);
    ylabel('Monostatic RCS (dBsm)', 'FontSize', 12);
    title(sprintf('Backscatter (Monostatic) vs Incident Angle  %s', compLabel), ...
          'FontSize', 11);
    grid on; box on;
    xlim([min(thetaInc) max(thetaInc)]);

    % Peak
    [peakVal, peakIdx] = max(backwardRCS);
    hold on;
    plot(thetaInc(peakIdx), peakVal, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    text(thetaInc(peakIdx), peakVal, ...
         sprintf('  %.1f dBsm @ %.0f°', peakVal, thetaInc(peakIdx)), ...
         'FontSize', 9, 'Color', 'r');

    % ---- Bistatic RCS at a few fixed incident angles ----
    subplot(1, 2, 2);
    hold on;
    N_show = min(5, N_inc);
    showIdx = unique(round(linspace(1, N_inc, N_show)));
    colors = lines(N_show);

    for i = 1:N_show
        idx = showIdx(i);
        rcsCut = rcsDB(:, idx);
        if size(rcsCut, 1) ~= length(thetaObs)
            rcsCut = rcsCut';
        end
        plot(thetaObs, rcsCut, '-', 'LineWidth', 1.8, 'Color', colors(i,:));

        % Mark incident direction on plot
        th_i = thetaInc(idx);
        yL = ylim;
        plot([th_i th_i], yL, ':', 'Color', colors(i,:), 'LineWidth', 1);
    end
    xlabel('Observation Angle \theta_o (deg)', 'FontSize', 12);
    ylabel('RCS (dBsm)', 'FontSize', 12);
    title(sprintf('Bistatic RCS at Selected Incident Angles  %s', compLabel), ...
          'FontSize', 11);
    grid on; box on;

    sgtitle(sprintf('Forward / Backward Scattering Analysis  —  %s  |  f = %.2f GHz', ...
            data.fileName, data.freq_ghz), 'FontSize', 12, 'Interpreter', 'none');

    if opt.saveFigs
        saveFig(fig, fullfile(opt.saveDir, [opt.figPrefix '_forward_backward']), opt);
    end
end

%% ========================================================================
%  FIGURE 5: 3D Surface of Bistatic RCS Matrix
% ========================================================================
function plotBistatic3D(thetaObs, thetaInc, rcsDB, compLabel, data, opt)
    fig = figure('Name', sprintf('FFE - 3D Bistatic [%s]', data.fileName), ...
                 'NumberTitle', 'off', ...
                 'Position', [50, 50, 900, 600], 'Color', 'w');

    [T_inc, T_obs] = meshgrid(thetaInc, thetaObs);

    % ---- 3D Surface ----
    subplot(1, 2, 1);
    surf(T_inc, T_obs, rcsDB, 'EdgeColor', 'interp', 'FaceColor', 'interp');
    xlabel('\theta_i (deg)', 'FontSize', 12);
    ylabel('\theta_o (deg)', 'FontSize', 12);
    zlabel('RCS (dBsm)', 'FontSize', 12);
    title(sprintf('Bistatic RCS Surface  %s', compLabel), 'FontSize', 11);
    colormap(gca, opt.colormap);
    colorbar;
    view([-45, 35]);
    if ~verLessThan('matlab', '8.4')
        lighting gouraud;
        camlight('headlight');
    end
    grid on;

    % ---- Contour ----
    subplot(1, 2, 2);
    nContours = 15;
    contourf(T_inc, T_obs, rcsDB, nContours, 'LineColor', 'none');
    hold on;
    contour(T_inc, T_obs, rcsDB, nContours, 'k-', 'LineWidth', 0.3);
    xlabel('\theta_i (deg)', 'FontSize', 12);
    ylabel('\theta_o (deg)', 'FontSize', 12);
    title(sprintf('Bistatic RCS Contour  %s', compLabel), 'FontSize', 11);
    colormap(gca, opt.colormap);
    colorbar;
    axis tight;

    sgtitle(sprintf('3D Bistatic RCS  —  %s  |  f = %.2f GHz  |  %s', ...
            data.fileName, data.freq_ghz, compLabel), ...
            'FontSize', 12, 'Interpreter', 'none');

    if opt.saveFigs
        saveFig(fig, fullfile(opt.saveDir, [opt.figPrefix '_bistatic_3d']), opt);
    end
end

%% ========================================================================
%  EXTRACTMONOSTATIC  Extract monostatic RCS from bistatic matrix
%  Finds the diagonal where theta_obs ≈ theta_inc
% ========================================================================
function [monoAngles, monoRCS] = extractMonostatic(thetaObs, thetaInc, rcsDB)
    % Find the overlapping angle range
    tMin = max(min(thetaObs), min(thetaInc));
    tMax = min(max(thetaObs), max(thetaInc));

    % Build a common angle axis
    N_mono = min(length(thetaObs), length(thetaInc));
    monoAngles = linspace(tMin, tMax, N_mono)';

    % Interpolate to get RCS at theta_obs = theta_inc = monoAngles
    monoRCS = zeros(N_mono, 1);
    for i = 1:N_mono
        th = monoAngles(i);
        % Find nearest observation and incident indices
        [~, obsIdx] = min(abs(thetaObs - th));
        [~, incIdx] = min(abs(thetaInc - th));
        if obsIdx <= size(rcsDB, 1) && incIdx <= size(rcsDB, 2)
            monoRCS(i) = rcsDB(obsIdx, incIdx);
        else
            monoRCS(i) = NaN;
        end
    end

    % Remove NaN
    valid = ~isnan(monoRCS);
    monoAngles = monoAngles(valid);
    monoRCS = monoRCS(valid);
end

%% ========================================================================
function saveFig(fig, basePath, opt)
    if ~exist(opt.saveDir, 'dir')
        [status, msg] = mkdir(opt.saveDir);
        if ~status
            warning('ffe_plot_bistatic:CannotCreateDir', ...
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
            warning('ffe_plot_bistatic:SaveFailed', ...
                    'Failed to save figure: %s', ME.message);
        end
    end
end
