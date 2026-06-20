% FFE_EXPORT  Export RCS data and figures to various formats
%   RCS数据与图形导出工具
%
%   Exports RCS data and/or figures to:
%     - PNG  (high-resolution raster, 300 DPI)
%     - SVG  (vector graphics)
%     - EPS  (encapsulated PostScript)
%     - FIG  (MATLAB figure format)
%     - CSV  (comma-separated values table)
%     - XLSX (Excel spreadsheet, multi-sheet)
%     - MAT  (MATLAB binary data)
%
%   Input:
%       data     - data struct from ffe_load_data (required for data export)
%       format   - string or cell array of strings:
%                  'png', 'svg', 'eps', 'fig', 'csv', 'xlsx', 'mat', 'all'
%       options  - (optional) struct with fields:
%           .saveDir    - export directory (default = '../results/ffe_export/')
%           .filePrefix - prefix for exported file names
%           .resolution - DPI for raster export (default = 300)
%           .exportData - true/false to export data files (default = true)
%           .exportFigs - true/false to export current figures (default = false)
%
%   Usage:
%     >> data = ffe_load_data('results/temp.dat');
%     >> ffe_export(data, {'png', 'csv'});
%     >> ffe_export(data, 'all', struct('saveDir', 'my_exports/'));
%
%   See also: ffe_main, ffe_load_data

function ffe_export(data, format, options)

    %% ---- Validate ----
    if nargin < 1 || isempty(data)
        error('ffe_export:NoData', 'Data struct required. Load with ffe_load_data first.');
    end
    if nargin < 2 || isempty(format)
        format = 'csv';
    end
    if nargin < 3, options = struct(); end

    opt = setDefaults(options, data);

    %% ---- Normalize format to cell array ----
    if ischar(format)
        if strcmpi(format, 'all')
            format = {'png', 'svg', 'csv', 'mat'};
        else
            format = {format};
        end
    end

    %% ---- Ensure export directory exists ----
    if ~exist(opt.saveDir, 'dir')
        [status, msg] = mkdir(opt.saveDir);
        if ~status
            error('ffe_export:CannotCreateDir', ...
                  'Cannot create export directory "%s": %s', opt.saveDir, msg);
        end
    end

    fprintf('========================================\n');
    fprintf('  FFE Export\n');
    fprintf('========================================\n');
    fprintf('  Source: %s\n', data.fileName);
    fprintf('  Target: %s\n', opt.saveDir);
    fprintf('  Formats: %s\n', strjoin(format, ', '));

    %% ---- Process each format ----
    for i = 1:length(format)
        fmt = lower(strtrim(format{i}));
        switch fmt
            case 'csv'
                exportCSV(data, opt);
            case 'xlsx'
                exportXLSX(data, opt);
            case 'mat'
                exportMAT(data, opt);
            case {'png', 'svg', 'eps', 'fig'}
                if opt.exportFigs
                    exportFigures(fmt, opt);
                else
                    fprintf('  Skipping %s (exportFigs=false)\n', fmt);
                end
            otherwise
                warning('ffe_export:UnknownFormat', ...
                        'Unknown export format: "%s". Skipping.', fmt);
        end
    end

    fprintf('========================================\n');
    fprintf('  Export complete.\n');
    fprintf('========================================\n');
end

%% ========================================================================
function opt = setDefaults(options, data)
    opt.saveDir    = getOpt(options, 'saveDir', ...
                             fullfile('..', 'results', 'ffe_export'));
    opt.filePrefix = getOpt(options, 'filePrefix', data.fileName);
    opt.resolution = getOpt(options, 'resolution', 300);
    opt.exportData = getOpt(options, 'exportData', true);
    opt.exportFigs = getOpt(options, 'exportFigs', false);
end

function val = getOpt(s, field, defaultVal)
    if isfield(s, field) && ~isempty(s.(field))
        val = s.(field);
    else
        val = defaultVal;
    end
end

%% ========================================================================
%  EXPORTCSV  Export RCS data as CSV file
% ========================================================================
function exportCSV(data, opt)
    if ~opt.exportData, return; end

    filePath = fullfile(opt.saveDir, [opt.filePrefix '.csv']);
    fprintf('  Exporting CSV: %s\n', filePath);

    fid = fopen(filePath, 'w');
    if fid == -1
        warning('ffe_export:CannotWrite', 'Cannot write to: %s', filePath);
        return;
    end

    % Header
    fprintf(fid, '# FFE Export — RCS Data\n');
    fprintf(fid, '# Source: %s\n', data.fileName);
    fprintf(fid, '# Type: %s\n', data.type);
    fprintf(fid, '# Grid: %d phi x %d theta\n', data.ip, data.it);
    fprintf(fid, '#\n');

    if data.is1D
        % 1D cut: single-column format
        fprintf(fid, 'Angle_deg,RCS_Theta_dBsm,RCS_Phi_dBsm\n');
        if data.ip == 1
            angle = data.theta(1, :);
        else
            angle = data.phi(:, 1)';
        end
        Sth = data.Sth(:)';
        Sph = data.Sph(:)';

        for i = 1:length(angle)
            fprintf(fid, '%.4f,%.6f,%.6f\n', angle(i), Sth(i), Sph(i));
        end
    else
        % 2D grid: theta as columns, phi as rows, separate blocks
        fprintf(fid, '# Theta (deg) — columns\n');
        fprintf(fid, 'Phi_deg\\Theta_deg');
        for j = 1:size(data.theta, 2)
            fprintf(fid, ',%.4f', data.theta(1, j));
        end
        fprintf(fid, '\n');

        fprintf(fid, '# RCS Theta (dBsm)\n');
        for i = 1:size(data.Sth, 1)
            fprintf(fid, '%.4f', data.phi(i, 1));
            for j = 1:size(data.Sth, 2)
                fprintf(fid, ',%.6f', data.Sth(i, j));
            end
            fprintf(fid, '\n');
        end

        fprintf(fid, '\n# RCS Phi (dBsm)\n');
        fprintf(fid, 'Phi_deg\\Theta_deg');
        for j = 1:size(data.theta, 2)
            fprintf(fid, ',%.4f', data.theta(1, j));
        end
        fprintf(fid, '\n');
        for i = 1:size(data.Sph, 1)
            fprintf(fid, '%.4f', data.phi(i, 1));
            for j = 1:size(data.Sph, 2)
                fprintf(fid, ',%.6f', data.Sph(i, j));
            end
            fprintf(fid, '\n');
        end
    end

    fclose(fid);
    fprintf('    Done: %s\n', filePath);
end

%% ========================================================================
%  EXPORTXLSX  Export RCS data as Excel spreadsheet
% ========================================================================
function exportXLSX(data, opt)
    if ~opt.exportData, return; end

    filePath = fullfile(opt.saveDir, [opt.filePrefix '.xlsx']);
    fprintf('  Exporting XLSX: %s\n', filePath);

    try
        if data.is1D
            % 1D cut: simple table
            if data.ip == 1
                angle = data.theta(1, :);
            else
                angle = data.phi(:, 1)';
            end
            T = table(angle(:), data.Sth(:), data.Sph(:), ...
                      'VariableNames', {'Angle_deg', 'RCS_Theta_dBsm', 'RCS_Phi_dBsm'});
            writetable(T, filePath, 'Sheet', 'RCS_Data');
        else
            % 2D grid: separate sheets for theta and phi components
            thetaRow = data.theta(1, :);
            T_th = array2table(data.Sth, ...
                   'VariableNames', arrayfun(@(x) sprintf('Th_%.1f', x), ...
                   thetaRow, 'UniformOutput', false));
            T_th = addvars(T_th, data.phi(:,1), 'Before', 1, ...
                   'NewVariableNames', {'Phi_deg'});
            writetable(T_th, filePath, 'Sheet', 'RCS_Theta');

            T_ph = array2table(data.Sph, ...
                   'VariableNames', arrayfun(@(x) sprintf('Th_%.1f', x), ...
                   thetaRow, 'UniformOutput', false));
            T_ph = addvars(T_ph, data.phi(:,1), 'Before', 1, ...
                   'NewVariableNames', {'Phi_deg'});
            writetable(T_ph, filePath, 'Sheet', 'RCS_Phi');
        end

        % Add metadata sheet
        metaTable = table({data.fileName; data.type; sprintf('%d×%d', data.ip, data.it)}, ...
                          'VariableNames', {'Value'}, ...
                          'RowNames', {'Source', 'Type', 'Grid'});
        writetable(metaTable, filePath, 'Sheet', 'Metadata', 'WriteRowNames', true);

        fprintf('    Done: %s\n', filePath);

    catch ME
        warning('ffe_export:XLSXFailed', ...
                'XLSX export failed: %s. Try installing/updating the MATLAB Excel support.', ...
                ME.message);
    end
end

%% ========================================================================
%  EXPORTMAT  Export data as MATLAB .mat file
% ========================================================================
function exportMAT(data, opt)
    if ~opt.exportData, return; end

    filePath = fullfile(opt.saveDir, [opt.filePrefix '_export.mat']);
    fprintf('  Exporting MAT: %s\n', filePath);

    % Extract all relevant fields for saving
    exportData = struct();
    exportData.type      = data.type;
    exportData.fileName  = data.fileName;
    exportData.theta     = data.theta;
    exportData.phi       = data.phi;
    exportData.Sth       = data.Sth;
    exportData.Sph       = data.Sph;
    exportData.ip        = data.ip;
    exportData.it        = data.it;
    exportData.N_angles  = data.N_angles;
    exportData.param     = data.param;

    if data.hasFreq
        exportData.S_complex = data.S_complex;
        exportData.freq      = data.freq;
        exportData.freq_ghz  = data.freq_ghz;
        exportData.N_f       = data.N_f;
    end

    save(filePath, '-struct', 'exportData', '-v7.3');
    fprintf('    Done: %s\n', filePath);
end

%% ========================================================================
%  EXPORTFIGURES  Export all open figure windows
% ========================================================================
function exportFigures(fmt, opt)
    figHandles = findobj('Type', 'figure');
    if isempty(figHandles)
        fprintf('  No figures to export.\n');
        return;
    end

    for i = 1:length(figHandles)
        fig = figHandles(i);
        figName = get(fig, 'Name');
        if isempty(figName)
            figName = sprintf('figure_%d', fig.Number);
        end
        % Sanitize filename
        safeName = matlab.lang.makeValidName(figName);

        filePath = fullfile(opt.saveDir, ...
                   sprintf('%s_fig%d_%s.%s', opt.filePrefix, fig.Number, safeName, fmt));
        fprintf('  Exporting figure: %s\n', filePath);

        switch fmt
            case 'png'
                try
                    exportgraphics(fig, filePath, 'Resolution', opt.resolution);
                catch
                    saveas(fig, filePath);
                end
            case 'svg'
                try
                    % SVG needs vector-friendly renderer
                    set(fig, 'Renderer', 'painters');
                    saveas(fig, filePath);
                catch ME
                    warning('ffe_export:SVGFailed', ...
                            'SVG export failed: %s', ME.message);
                end
            case 'eps'
                try
                    set(fig, 'Renderer', 'painters');
                    saveas(fig, filePath, 'epsc');
                catch ME
                    warning('ffe_export:EPSFailed', ...
                            'EPS export failed: %s', ME.message);
                end
            case 'fig'
                savefig(fig, filePath);
        end
    end

    fprintf('    Exported %d figure(s).\n', length(figHandles));
end
