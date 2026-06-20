% FFE_PLOT_MONOSTATIC  Monostatic RCS plotting from .ffe bistatic data
%   单站RCS绘图 — 从双站数据矩阵对角线提取
%
%   Extracts the monostatic RCS (theta_obs = theta_inc diagonal) from
%   POSTFEKO .ffe bistatic data and creates publication-quality plots.
%
%   Features:
%     - Cartesian plot with peak, -3dB beamwidth, and sidelobe markers
%     - Polar plot
%     - Multi-file comparison (overlay up to 4 datasets)
%     - Key metrics table (peak RCS, beamwidth, mean RCS, etc.)
%
%   Input:
%       filePaths - string or cell array of .ffe file paths
%                   If empty, opens file dialog (multi-select enabled)
%       options   - (optional) struct:
%           .component   - 'total' (default), 'theta', 'phi'
%           .normalize   - normalize to 0dB peak for comparison (default=false)
%           .saveFigs    - save figures to file (default=false)
%           .saveDir     - output directory (default='./')
%
%   Usage:
%     >> ffe_plot_monostatic('ffe_file/result_rocket_po.ffe');
%     >> ffe_plot_monostatic({'ffe_file/result_rocket_po.ffe', ...
%                              'ffe_file/result_rocket_lpo.ffe'});
%     >> ffe_plot_monostatic();  % opens file dialog
%
%   See also: ffe_load_data, ffe_plot_bistatic

function ffe_plot_monostatic(filePaths, options)

    %% ---- Handle input ----
    if nargin < 1 || isempty(filePaths)
        % Open file dialog with multi-select
        [fileNames, pathName] = uigetfile( ...
            {'*.ffe', 'POSTFEKO Far-Field Files (*.ffe)'}, ...
            'Select FFE file(s) for Monostatic Plot', ...
            'ffe_file/', 'MultiSelect', 'on');
        if isequal(fileNames, 0)
            disp('No file selected. Exiting.');
            return;
        end
        if ischar(fileNames)
            fileNames = {fileNames};
        end
        filePaths = cellfun(@(f) fullfile(pathName, f), fileNames, 'UniformOutput', false);
    end

    if ischar(filePaths) || isstring(filePaths)
        filePaths = {char(filePaths)};
    end

    if nargin < 2, options = struct(); end
    opt = setDefaults(options);

    N_files = length(filePaths);
    fprintf('=== Monostatic RCS Plot ===\n');
    fprintf('  %d file(s) to process\n', N_files);

    %% ---- Load all files and extract monostatic cuts ----
    monoData = cell(N_files, 1);
    fileNames = cell(N_files, 1);
    legendNames = cell(N_files, 1);

    for i = 1:N_files
        fprintf('  Loading [%d/%d]: %s\n', i, N_files, filePaths{i});
        data = ffe_load_data(filePaths{i});

        % Choose RCS component
        switch lower(opt.component)
            case 'theta'
                rcsLin = data.RCS_theta;
                compStr = 'Theta';
            case 'phi'
                rcsLin = data.RCS_phi;
                compStr = 'Phi';
            case 'total'
                rcsLin = data.RCS_total;
                compStr = 'Total';
        end

        rcsDB = 10 * log10(max(rcsLin, 1e-30));
        thetaObs = data.theta(1, :)';
        thetaInc = data.theta_inc(:);

        % Extract monostatic diagonal
        [angles, rcsMono] = extractMonostatic(thetaObs, thetaInc, rcsDB);

        monoData{i}.angles = angles;
        monoData{i}.rcs    = rcsMono;
        monoData{i}.freqGHz = data.freq_ghz;
        monoData{i}.name   = data.fileName;
        monoData{i}.source = data.source;

        [~, fname] = fileparts(filePaths{i});
        fileNames{i} = fname;
        legendNames{i} = sprintf('%s  (%.2f GHz)', fname, data.freq_ghz);
    end

    %% ---- Figure 1: Cartesian plot ----
    plotMonostaticCartesian(monoData, legendNames, fileNames, opt);

    %% ---- Figure 2: Polar plot ----
    plotMonostaticPolar(monoData, legendNames, fileNames, opt);

    %% ---- Print metrics table ----
    printMetrics(monoData, legendNames);
end

%% ========================================================================
function opt = setDefaults(options)
    opt.component = getOpt(options, 'component', 'total');
    opt.normalize = getOpt(options, 'normalize', false);
    opt.saveFigs  = getOpt(options, 'saveFigs', false);
    opt.saveDir   = getOpt(options, 'saveDir', './');
end

function val = getOpt(s, field, defVal)
    if isfield(s, field) && ~isempty(s.(field))
        val = s.(field);
    else
        val = defVal;
    end
end

%% ========================================================================
%  Cartesian monostatic RCS plot
%% ========================================================================
function plotMonostaticCartesian(monoData, legendNames, fileNames, opt)
    fig = figure('Name', 'Monostatic RCS - Cartesian', ...
                 'NumberTitle', 'off', ...
                 'Position', [50, 80, 900, 550], 'Color', 'w');

    N = length(monoData);
    colors = lines(N);

    hold on;
    h = zeros(N, 1);

    for i = 1:N
        ang = monoData{i}.angles;
        rcs = monoData{i}.rcs;

        if opt.normalize
            rcs = rcs - max(rcs);
            yLab = 'Normalized Monostatic RCS (dB)';
        else
            yLab = 'Monostatic RCS (dBsm)';
        end

        h(i) = plot(ang, rcs, '-', 'LineWidth', 2.2, 'Color', colors(i,:));

        % Peak marker
        [peakVal, peakIdx] = max(rcs);
        plot(ang(peakIdx), peakVal, 'v', 'MarkerSize', 10, ...
             'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k');
        text(ang(peakIdx), peakVal, ...
             sprintf(' %.0f°', ang(peakIdx)), ...
             'FontSize', 8, 'Color', colors(i,:), 'VerticalAlignment', 'bottom');

        % -3dB beamwidth (show only for first file to avoid clutter)
        if i == 1 && N == 1
            [bw, li, ri] = beamwidth(rcs, ang, 3);
            if ~isnan(bw)
                level3dB = peakVal - 3;
                yRange = ylim;
                plot([ang(1) ang(end)], [level3dB level3dB], '--', ...
                     'Color', [0.5 0.5 0.5], 'LineWidth', 1);
                text(mean([ang(1) ang(end)]), level3dB + 0.5, ...
                     sprintf('-3 dB BW: %.1f°', bw), ...
                     'FontSize', 10, 'Color', [0.4 0.4 0.4], ...
                     'HorizontalAlignment', 'center');
            end
        end
    end

    xlabel('\theta (deg)', 'FontSize', 13);
    ylabel(yLab, 'FontSize', 13);
    title('Monostatic RCS  (\theta_o = \theta_i)', 'FontSize', 13);
    legend(h, legendNames, 'Location', 'best', 'FontSize', 10, 'Interpreter', 'none');
    grid on; box on;
    xlim([min(monoData{1}.angles) max(monoData{1}.angles)]);

    if opt.saveFigs
        saveFig(fig, fullfile(opt.saveDir, 'monostatic_cartesian'), opt);
    end
end

%% ========================================================================
%  Polar monostatic RCS plot
%% ========================================================================
function plotMonostaticPolar(monoData, legendNames, fileNames, opt)
    fig = figure('Name', 'Monostatic RCS - Polar', ...
                 'NumberTitle', 'off', ...
                 'Position', [100, 50, 620, 580], 'Color', 'w');

    N = length(monoData);
    colors = lines(N);

    h = zeros(N, 1);

    for i = 1:N
        ang = monoData{i}.angles;
        rcs = monoData{i}.rcs;

        % Shift to positive range for polar
        if opt.normalize
            rcs = rcs - max(rcs);
        end
        minR = min(rcs);
        rPolar = rcs - minR + 5;

        if i == 1
            h(i) = polarplot(deg2rad(ang), rPolar, '-', 'LineWidth', 2.2, 'Color', colors(i,:));
        else
            hold on;
            h(i) = polarplot(deg2rad(ang), rPolar, '-', 'LineWidth', 2.2, 'Color', colors(i,:));
        end
    end

    title(sprintf('Monostatic RCS Polar  (\theta_o = \\theta_i)'), 'FontSize', 13);
    legend(h, legendNames, 'Location', 'bestoutside', 'FontSize', 9, 'Interpreter', 'none');

    if opt.saveFigs
        saveFig(fig, fullfile(opt.saveDir, 'monostatic_polar'), opt);
    end
end

%% ========================================================================
%  Print key metrics
%% ========================================================================
function printMetrics(monoData, legendNames)
    fprintf('\n========================================\n');
    fprintf('  Monostatic RCS Metrics\n');
    fprintf('========================================\n');
    fprintf('  %-30s %10s %10s %10s %10s\n', ...
            'File', 'Peak(dBsm)', 'Peak@deg', '-3dB BW°', 'Mean(dBsm)');
    fprintf('  %-30s %10s %10s %10s %10s\n', ...
            '----', '----------', '--------', '--------', '---------');

    for i = 1:length(monoData)
        ang = monoData{i}.angles;
        rcs = monoData{i}.rcs;

        [peakVal, peakIdx] = max(rcs);
        peakAng = ang(peakIdx);
        [bw, ~, ~] = beamwidth(rcs, ang, 3);
        meanRCS = mean(rcs);

        if isnan(bw)
            bwStr = 'N/A';
        else
            bwStr = sprintf('%.1f', bw);
        end

        fprintf('  %-30s %10.1f %10.1f %10s %10.1f\n', ...
                legendNames{i}, peakVal, peakAng, bwStr, meanRCS);
    end
    fprintf('========================================\n');
end

%% ========================================================================
%  Extract monostatic RCS (theta_obs = theta_inc diagonal)
%% ========================================================================
function [monoAngles, monoRCS] = extractMonostatic(thetaObs, thetaInc, rcsDB)
    tMin = max(min(thetaObs), min(thetaInc));
    tMax = min(max(thetaObs), max(thetaInc));
    N_mono = min(length(thetaObs), length(thetaInc));
    monoAngles = linspace(tMin, tMax, N_mono)';

    monoRCS = zeros(N_mono, 1);
    for i = 1:N_mono
        th = monoAngles(i);
        [~, obsIdx] = min(abs(thetaObs - th));
        [~, incIdx] = min(abs(thetaInc - th));
        if obsIdx <= size(rcsDB, 1) && incIdx <= size(rcsDB, 2)
            monoRCS(i) = rcsDB(obsIdx, incIdx);
        else
            monoRCS(i) = NaN;
        end
    end
    valid = ~isnan(monoRCS);
    monoAngles = monoAngles(valid);
    monoRCS = monoRCS(valid);
end

%% ========================================================================
function saveFig(fig, basePath, opt)
    if ~exist(opt.saveDir, 'dir')
        [status, msg] = mkdir(opt.saveDir);
        if ~status
            warning('Cannot create: %s', msg);
            return;
        end
    end
    try
        exportgraphics(fig, [basePath '.png'], 'Resolution', 300);
        fprintf('  Saved: %s.png\n', basePath);
    catch
        try
            saveas(fig, [basePath '.png']);
        catch ME
            warning('Save failed: %s', ME.message);
        end
    end
end
