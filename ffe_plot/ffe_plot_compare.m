% FFE_PLOT_COMPARE  Comparative RCS analysis between two datasets
%   双数据集RCS方向图对比分析
%
%   Loads two RCS result files and creates comparative visualizations:
%     - Overlay plot (both patterns on same axes)
%     - Difference plot (ΔdBsm vs angle)
%     - Side-by-side comparison with synchronized axes
%     - Statistical metrics: correlation coefficient, RMS difference,
%       peak difference, integrated RCS difference
%
%   Input:
%       data1   - data struct from ffe_load_data (reference/baseline)
%       data2   - data struct from ffe_load_data (comparison)
%       options - (optional) struct with fields:
%           .component    - 'theta' (default), 'phi', or 'both'
%           .normalize    - true to align peaks before comparison (default = false)
%           .saveFigs     - true/false to save figures (default = false)
%           .saveDir      - directory for saved figures
%           .figPrefix    - prefix for saved figure filenames
%
%   Usage:
%     >> d1 = ffe_load_data('results/run1.dat');
%     >> d2 = ffe_load_data('results/run2.dat');
%     >> ffe_plot_compare(d1, d2);
%     >> ffe_plot_compare(d1, d2, struct('normalize', true, 'component', 'both'));
%
%   See also: ffe_main, ffe_load_data, ffe_plot_rcs_2d

function ffe_plot_compare(data1, data2, options)

    %% ---- Validate ----
    if nargin < 2
        error('ffe_plot_compare:NeedTwoDatasets', ...
              'Two datasets required. Load them with ffe_load_data first.');
    end
    if nargin < 3, options = struct(); end
    opt = setDefaults(options);

    %% ---- Extract 1D slices for comparison ----
    [ang1, rTh1, rPh1, angLab1, sliceLab1] = extractSlice(data1);
    [ang2, rTh2, rPh2, angLab2, sliceLab2] = extractSlice(data2);

    % For comparison, both datasets need compatible angle grids
    % Use the finer grid and interpolate the coarser one
    [angle, rTh1_i, rTh2_i, rPh1_i, rPh2_i, angLabel] = ...
        alignAngleGrids(ang1, rTh1, rPh1, ang2, rTh2, rPh2, angLab1, angLab2);

    %% ---- Normalize if requested ----
    if opt.normalize
        rTh1_i = rTh1_i - max(rTh1_i);
        rTh2_i = rTh2_i - max(rTh2_i);
        rPh1_i = rPh1_i - max(rPh1_i);
        rPh2_i = rPh2_i - max(rPh2_i);
        yLabelStr = 'Normalized RCS (dB)';
    else
        yLabelStr = 'RCS (dBsm)';
    end

    %% ---- Figure 1: Overlay comparison ----
    plotOverlay(angle, rTh1_i, rTh2_i, rPh1_i, rPh2_i, angLabel, ...
                yLabelStr, data1, data2, opt);

    %% ---- Figure 2: Difference plot ----
    plotDifference(angle, rTh1_i, rTh2_i, rPh1_i, rPh2_i, angLabel, ...
                   data1, data2, opt);

    %% ---- Figure 3: Side-by-side comparison ----
    plotSideBySide(angle, rTh1_i, rTh2_i, rPh1_i, rPh2_i, angLabel, ...
                   yLabelStr, data1, data2, opt);

    %% ---- Print statistics ----
    printStatistics(angle, rTh1_i, rTh2_i, rPh1_i, rPh2_i, data1, data2, opt);
end

%% ========================================================================
function opt = setDefaults(options)
    opt.component = getOpt(options, 'component', 'theta');
    opt.normalize = getOpt(options, 'normalize', false);
    opt.saveFigs  = getOpt(options, 'saveFigs', false);
    opt.saveDir   = getOpt(options, 'saveDir', fullfile('..', 'results'));
    opt.figPrefix = getOpt(options, 'figPrefix', 'comparison');
end

function val = getOpt(s, field, defaultVal)
    if isfield(s, field) && ~isempty(s.(field))
        val = s.(field);
    else
        val = defaultVal;
    end
end

%% ========================================================================
function [angle, rTh, rPh, angLab, sliceLab] = extractSlice(data)
    if data.ip == 1
        angle = data.theta(1, :);
        rTh   = data.Sth(1, :);
        rPh   = data.Sph(1, :);
        angLab   = '\theta (deg)';
        sliceLab = sprintf('\\phi=%.0f°', data.phi(1,1));
    elseif data.it == 1
        angle = data.phi(:, 1)';
        rTh   = data.Sth(:, 1)';
        rPh   = data.Sph(:, 1)';
        angLab   = '\phi (deg)';
        sliceLab = sprintf('\\theta=%.0f°', data.theta(1,1));
    else
        midPhi = round(data.ip / 2);
        angle = data.theta(midPhi, :);
        rTh   = data.Sth(midPhi, :);
        rPh   = data.Sph(midPhi, :);
        angLab   = '\theta (deg)';
        sliceLab = sprintf('\\phi=%.0f° (center)', data.phi(midPhi,1));
    end
    angle = angle(:)';
    rTh = rTh(:)';
    rPh = rPh(:)';
end

%% ========================================================================
function [angle, rTh1, rTh2, rPh1, rPh2, angLabel] = ...
        alignAngleGrids(ang1, rTh1, rPh1, ang2, rTh2, rPh2, angLab1, ~)

    angLabel = angLab1;

    % Use union of angle ranges
    aMin = max(min(ang1), min(ang2));
    aMax = min(max(ang1), max(ang2));
    N = max(length(ang1), length(ang2));
    angle = linspace(aMin, aMax, N);

    rTh1 = interp1(ang1, rTh1, angle, 'linear', 'extrap');
    rPh1 = interp1(ang1, rPh1, angle, 'linear', 'extrap');
    rTh2 = interp1(ang2, rTh2, angle, 'linear', 'extrap');
    rPh2 = interp1(ang2, rPh2, angle, 'linear', 'extrap');
end

%% ========================================================================
%  FIGURE 1: Overlay plot
% ========================================================================
function plotOverlay(angle, rTh1, rTh2, rPh1, rPh2, angLabel, ...
                     yLabel, data1, data2, opt)
    fig = figure('Name', 'FFE - Comparison Overlay', 'NumberTitle', 'off', ...
                 'Position', [60, 80, 950, 520], 'Color', 'w');

    colors = lines(4);

    subplot(1, 2, 1);
    hold on;
    h1 = plot(angle, rTh1, '-', 'LineWidth', 2.2, 'Color', colors(1,:));
    h2 = plot(angle, rTh2, '--', 'LineWidth', 2.2, 'Color', colors(2,:));
    xlabel(angLabel, 'FontSize', 12);
    ylabel(yLabel, 'FontSize', 12);
    title('S_{\theta} Comparison', 'FontSize', 12);
    legend([h1, h2], {data1.fileName, data2.fileName}, ...
           'Location', 'best', 'FontSize', 10);
    grid on; box on;
    xlim([min(angle) max(angle)]);

    subplot(1, 2, 2);
    hold on;
    h3 = plot(angle, rPh1, '-', 'LineWidth', 2.2, 'Color', colors(3,:));
    h4 = plot(angle, rPh2, '--', 'LineWidth', 2.2, 'Color', colors(4,:));
    xlabel(angLabel, 'FontSize', 12);
    ylabel(yLabel, 'FontSize', 12);
    title('S_{\phi} Comparison', 'FontSize', 12);
    legend([h3, h4], {data1.fileName, data2.fileName}, ...
           'Location', 'best', 'FontSize', 10);
    grid on; box on;
    xlim([min(angle) max(angle)]);

    sgtitle(sprintf('RCS Comparison Overlay  —  %s  vs  %s', ...
            data1.fileName, data2.fileName), 'FontSize', 12, 'Interpreter', 'none');

    if opt.saveFigs
        saveFig(fig, fullfile(opt.saveDir, [opt.figPrefix '_overlay']), opt);
    end
end

%% ========================================================================
%  FIGURE 2: Difference plot
% ========================================================================
function plotDifference(angle, rTh1, rTh2, rPh1, rPh2, angLabel, ...
                        data1, data2, opt)
    fig = figure('Name', 'FFE - Difference Plot', 'NumberTitle', 'off', ...
                 'Position', [60, 60, 900, 500], 'Color', 'w');

    diffTh = rTh2 - rTh1;
    diffPh = rPh2 - rPh1;

    subplot(1, 2, 1);
    plot(angle, diffTh, 'b-', 'LineWidth', 2);
    hold on;
    plot([min(angle) max(angle)], [0 0], 'k--', 'LineWidth', 0.8);

    % Highlight positive/negative regions
    posIdx = diffTh >= 0;
    if any(posIdx)
        area(angle(posIdx), diffTh(posIdx), 'FaceColor', 'r', ...
             'FaceAlpha', 0.15, 'EdgeColor', 'none');
    end
    negIdx = diffTh < 0;
    if any(negIdx)
        area(angle(negIdx), diffTh(negIdx), 'FaceColor', 'b', ...
             'FaceAlpha', 0.15, 'EdgeColor', 'none');
    end

    xlabel(angLabel, 'FontSize', 12);
    ylabel('\Delta S_{\theta} (dB)', 'FontSize', 12);
    title(sprintf(['\\DeltaS_{\\theta}  (%s - %s)'], data2.fileName, data1.fileName), ...
          'FontSize', 11, 'Interpreter', 'none');
    grid on; box on;
    xlim([min(angle) max(angle)]);

    subplot(1, 2, 2);
    plot(angle, diffPh, 'r-', 'LineWidth', 2);
    hold on;
    plot([min(angle) max(angle)], [0 0], 'k--', 'LineWidth', 0.8);

    posIdx = diffPh >= 0;
    if any(posIdx)
        area(angle(posIdx), diffPh(posIdx), 'FaceColor', 'r', ...
             'FaceAlpha', 0.15, 'EdgeColor', 'none');
    end
    negIdx = diffPh < 0;
    if any(negIdx)
        area(angle(negIdx), diffPh(negIdx), 'FaceColor', 'b', ...
             'FaceAlpha', 0.15, 'EdgeColor', 'none');
    end

    xlabel(angLabel, 'FontSize', 12);
    ylabel('\Delta S_{\phi} (dB)', 'FontSize', 12);
    title(sprintf(['\\DeltaS_{\\phi}  (%s - %s)'], data2.fileName, data1.fileName), ...
          'FontSize', 11, 'Interpreter', 'none');
    grid on; box on;
    xlim([min(angle) max(angle)]);

    % RMS difference annotation
    rmsTh = sqrt(mean(diffTh.^2));
    rmsPh = sqrt(mean(diffPh.^2));
    annotation('textbox', [0.35 0.02 0.3 0.06], ...
               'String', sprintf('RMS ΔS_θ: %.2f dB  |  RMS ΔS_ϕ: %.2f dB', rmsTh, rmsPh), ...
               'FontSize', 10, 'HorizontalAlignment', 'center', ...
               'BackgroundColor', 'w', 'EdgeColor', [0.5 0.5 0.5]);

    sgtitle(sprintf('RCS Difference  —  %s  vs  %s', ...
            data1.fileName, data2.fileName), 'FontSize', 12, 'Interpreter', 'none');

    if opt.saveFigs
        saveFig(fig, fullfile(opt.saveDir, [opt.figPrefix '_difference']), opt);
    end
end

%% ========================================================================
%  FIGURE 3: Side-by-side comparison
% ========================================================================
function plotSideBySide(angle, rTh1, rTh2, rPh1, rPh2, angLabel, ...
                        yLabel, data1, data2, opt)
    fig = figure('Name', 'FFE - Side by Side', 'NumberTitle', 'off', ...
                 'Position', [40, 50, 1050, 420], 'Color', 'w');

    colors = lines(4);

    % ---- Dataset 1 ----
    subplot(1, 2, 1);
    hold on;
    plot(angle, rTh1, '-', 'LineWidth', 2, 'Color', colors(1,:));
    plot(angle, rPh1, '--', 'LineWidth', 2, 'Color', colors(2,:));
    xlabel(angLabel, 'FontSize', 12);
    ylabel(yLabel, 'FontSize', 12);
    title(sprintf('Dataset 1: %s', data1.fileName), 'FontSize', 11, 'Interpreter', 'none');
    legend('S_{\theta}', 'S_{\phi}', 'Location', 'best', 'FontSize', 10);
    grid on; box on;
    xlim([min(angle) max(angle)]);

    % ---- Dataset 2 ----
    subplot(1, 2, 2);
    hold on;
    plot(angle, rTh2, '-', 'LineWidth', 2, 'Color', colors(3,:));
    plot(angle, rPh2, '--', 'LineWidth', 2, 'Color', colors(4,:));
    xlabel(angLabel, 'FontSize', 12);
    ylabel(yLabel, 'FontSize', 12);
    title(sprintf('Dataset 2: %s', data2.fileName), 'FontSize', 11, 'Interpreter', 'none');
    legend('S_{\theta}', 'S_{\phi}', 'Location', 'best', 'FontSize', 10);
    grid on; box on;
    xlim([min(angle) max(angle)]);

    sgtitle('Side-by-Side RCS Comparison', 'FontSize', 12);

    if opt.saveFigs
        saveFig(fig, fullfile(opt.saveDir, [opt.figPrefix '_sidebyside']), opt);
    end
end

%% ========================================================================
%  STATISTICS  Print quantitative comparison metrics
% ========================================================================
function printStatistics(angle, rTh1, rTh2, rPh1, rPh2, data1, data2, opt)
    fprintf('\n========================================\n');
    fprintf('  Comparison Statistics\n');
    fprintf('========================================\n');
    fprintf('  Reference:  %s\n', data1.fileName);
    fprintf('  Comparison: %s\n', data2.fileName);
    fprintf('----------------------------------------\n');

    % Correlation coefficients
    ccTh = corrcoef(rTh1(:), rTh2(:));
    ccTh = ccTh(1, 2);
    ccPh = corrcoef(rPh1(:), rPh2(:));
    ccPh = ccPh(1, 2);

    % RMS difference
    rmsTh = sqrt(mean((rTh2 - rTh1).^2));
    rmsPh = sqrt(mean((rPh2 - rPh1).^2));

    % Peak difference
    peakDiffTh = max(rTh2) - max(rTh1);
    peakDiffPh = max(rPh2) - max(rPh1);

    % Mean difference (bias)
    meanDiffTh = mean(rTh2 - rTh1);
    meanDiffPh = mean(rPh2 - rPh1);

    % Integrated RCS difference
    rcs1Th_lin = 10.^(rTh1 / 10);
    rcs2Th_lin = 10.^(rTh2 / 10);
    intDiffTh = 10 * log10(sum(rcs2Th_lin) / max(sum(rcs1Th_lin), 1e-30));

    rcs1Ph_lin = 10.^(rPh1 / 10);
    rcs2Ph_lin = 10.^(rPh2 / 10);
    intDiffPh = 10 * log10(sum(rcs2Ph_lin) / max(sum(rcs1Ph_lin), 1e-30));

    fprintf('  S_theta:\n');
    fprintf('    Correlation Coeff:     %+.4f\n', ccTh);
    fprintf('    RMS Difference:        %.3f dB\n', rmsTh);
    fprintf('    Mean Bias:             %+.3f dB\n', meanDiffTh);
    fprintf('    Peak Difference:       %+.3f dB\n', peakDiffTh);
    fprintf('    Integrated RCS Diff:   %+.3f dB\n', intDiffTh);
    fprintf('  S_phi:\n');
    fprintf('    Correlation Coeff:     %+.4f\n', ccPh);
    fprintf('    RMS Difference:        %.3f dB\n', rmsPh);
    fprintf('    Mean Bias:             %+.3f dB\n', meanDiffPh);
    fprintf('    Peak Difference:       %+.3f dB\n', peakDiffPh);
    fprintf('    Integrated RCS Diff:   %+.3f dB\n', intDiffPh);
    fprintf('========================================\n');
end

%% ========================================================================
function saveFig(fig, basePath, opt)
    if ~exist(opt.saveDir, 'dir')
        [status, msg] = mkdir(opt.saveDir);
        if ~status
            warning('ffe_plot_compare:CannotCreateDir', ...
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
            warning('ffe_plot_compare:SaveFailed', ...
                    'Failed to save figure: %s', ME.message);
        end
    end
end
