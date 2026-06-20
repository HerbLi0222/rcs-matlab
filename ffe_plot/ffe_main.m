% FFE_MAIN  FFE Post-Processing Toolbox — Interactive Main Menu
%   FFE（远场电磁）后处理工具箱 — 交互式主界面
%
%   A unified post-processing interface for RCS simulation results.
%   Supports POSTFEKO .ffe, .dat (monostatic/bistatic), and .mat (wideband) files.
%
%   Features:
%     1. Load Data          — Load .ffe / .dat / .mat RCS result files
%     2. 2D RCS Plot        — Enhanced cartesian + polar RCS plots
%     3. 3D RCS Plot        — Spherical surface, mesh, heatmap visualization
%     4. Frequency Analysis — Mag/phase, quick range profile, spectrogram
%     5. Bistatic Analysis  — Heatmap, monostatic cut, incident-angle cuts
%     6. Compare Datasets   — Overlay, difference, side-by-side, statistics
%     7. Export             — Export to CSV, XLSX, PNG, SVG, MAT
%     8. Batch Process      — Process all files in a directory
%     0. Exit
%
%   Usage:
%     >> ffe_main                               % Interactive menu mode
%     >> ffe_main('ffe_file/result_rocket.ffe')  % Load a specific file
%     >> ffe_main('ffe_file', 'batch')           % Batch process a directory
%
%   See also: ffe_load_data, ffe_plot_rcs_2d, ffe_plot_rcs_3d,
%             ffe_plot_bistatic, ffe_plot_frequency, ffe_plot_compare, ffe_export

function ffe_main(varargin)

    %% ---- Setup paths ----
    ffeDir = fileparts(mfilename('fullpath'));
    addpath(ffeDir);
    libDir = fullfile(ffeDir, '..', 'lib');
    if exist(libDir, 'dir')
        addpath(libDir);
    end

    %% ---- Parse input arguments ----
    if nargin >= 1 && ~isempty(varargin{1})
        arg1 = varargin{1};
        if isfolder(arg1)
            batchProcess(arg1);
            return;
        elseif isfile(arg1)
            try
                data = ffe_load_data(arg1);
                runInteractiveMenu(data);
            catch ME
                fprintf(2, 'Error loading file: %s\n', ME.message);
            end
            return;
        end
    end

    %% ---- Interactive menu mode ----
    data = [];
    runInteractiveMenu(data);
end

%% ========================================================================
function runInteractiveMenu(data)
    while true
        printHeader();
        printMenu(data);

        choice = input('  Enter choice (0-8): ', 's');
        choice = strtrim(choice);

        switch choice
            case '1'
                % ---- Load Data ----
                try
                    data = ffe_load_data();
                catch ME
                    fprintf(2, '  Error: %s\n', ME.message);
                end

            case '2'
                % ---- 2D RCS Plot ----
                if ~checkData(data), continue; end
                try
                    if isfield(data, 'bistatic') && data.bistatic
                        data2d = selectBistaticBlock(data, '2D plot');
                        ffe_plot_rcs_2d(data2d, getPlotOptions());
                    else
                        ffe_plot_rcs_2d(data, getPlotOptions());
                    end
                catch ME
                    fprintf(2, '  Error: %s\n', ME.message);
                end

            case '3'
                % ---- 3D RCS Plot ----
                if ~checkData(data), continue; end
                try
                    if isfield(data, 'bistatic') && data.bistatic
                        data3d = selectBistaticBlock(data, '3D plot');
                        ffe_plot_rcs_3d(data3d, get3DOptions());
                    else
                        if ~data.is2D
                            fprintf('\n  Warning: Data is 1D. Only heatmap will be shown.\n');
                        end
                        ffe_plot_rcs_3d(data, get3DOptions());
                    end
                catch ME
                    fprintf(2, '  Error: %s\n', ME.message);
                end

            case '4'
                % ---- Frequency Analysis ----
                if ~checkData(data), continue; end
                if ~data.hasFreq
                    fprintf('\n  Error: No frequency-domain data. Load a wideband_scattering_*.mat file.\n');
                    continue;
                end
                try
                    ffe_plot_frequency(data, getFreqOptions(data));
                catch ME
                    fprintf(2, '  Error: %s\n', ME.message);
                end

            case '5'
                % ---- Bistatic Analysis (NEW) ----
                if ~checkData(data), continue; end
                if ~isfield(data, 'bistatic') || ~data.bistatic
                    fprintf('\n  Error: Not bistatic data. Load a .ffe file first.\n');
                    continue;
                end
                try
                    ffe_plot_bistatic(data, getBistaticOptions(data));
                catch ME
                    fprintf(2, '  Error: %s\n', ME.message);
                end

            case '6'
                % ---- Compare Datasets ----
                fprintf('\n  --- Load Reference Dataset ---\n');
                try
                    d1 = ffe_load_data();
                catch ME
                    fprintf(2, '  Error: %s\n', ME.message);
                    continue;
                end
                fprintf('\n  --- Load Comparison Dataset ---\n');
                try
                    d2 = ffe_load_data();
                catch ME
                    fprintf(2, '  Error: %s\n', ME.message);
                    continue;
                end
                try
                    if isfield(d1, 'bistatic') && d1.bistatic
                        d1s = selectBistaticBlock(d1, 'comparison (reference)');
                    else
                        d1s = d1;
                    end
                    if isfield(d2, 'bistatic') && d2.bistatic
                        d2s = selectBistaticBlock(d2, 'comparison (target)');
                    else
                        d2s = d2;
                    end
                    ffe_plot_compare(d1s, d2s, getCompareOptions());
                catch ME
                    fprintf(2, '  Error: %s\n', ME.message);
                end

            case '7'
                % ---- Export ----
                if ~checkData(data), continue; end
                try
                    fmt = getExportFormat();
                    ffe_export(data, fmt, getExportOptions());
                catch ME
                    fprintf(2, '  Error: %s\n', ME.message);
                end

            case '8'
                % ---- Batch Process ----
                fprintf('  Default: ./ffe_file/\n');
                dirPath = input('  Enter directory path: ', 's');
                dirPath = strtrim(dirPath);
                if isempty(dirPath)
                    dirPath = fullfile(fileparts(mfilename('fullpath')), 'ffe_file');
                end
                if isfolder(dirPath)
                    batchProcess(dirPath);
                else
                    fprintf(2, '  Directory not found: %s\n', dirPath);
                end

            case '0'
                % ---- Exit ----
                fprintf('\n  Exiting FFE Post-Processing Toolbox. Goodbye!\n\n');
                break;

            otherwise
                fprintf(2, '  Invalid choice. Please enter 0-8.\n');
        end

        if ~isempty(choice) && choice ~= '0'
            fprintf('\n');
            input('  Press Enter to continue...', 's');
        end
    end
end

%% ========================================================================
function printHeader()
    clc;
    fprintf('========================================\n');
    fprintf('    FFE Post-Processing Toolbox v1.1\n');
    fprintf('    远场电磁后处理工具箱\n');
    fprintf('    Supports: POSTFEKO .ffe / .dat / .mat\n');
    fprintf('========================================\n');
end

function printMenu(data)
    if isempty(data)
        statusStr = '[No data loaded]';
    else
        if isfield(data, 'bistatic') && data.bistatic
            statusStr = sprintf('[Loaded: %s | %s Bistatic | %d obs × %d inc | f=%.2f GHz]', ...
                data.fileName, data.type, data.N_angles, data.N_inc, data.freq_ghz);
        else
            statusStr = sprintf('[Loaded: %s | %s | %d×%d grid | %d angles]', ...
                data.fileName, data.type, data.ip, data.it, data.N_angles);
            if data.hasFreq
                statusStr = [statusStr sprintf(' | %d freq points', data.N_f)]; %#ok<AGROW>
            end
        end
    end
    fprintf('  %s\n\n', statusStr);

    fprintf('  --- Analysis ---\n');
    fprintf('  1. Load Data (.ffe / .dat / .mat)\n');
    fprintf('  2. 2D RCS Plot (Cartesian + Polar)\n');
    fprintf('  3. 3D RCS Plot (Spherical / Mesh / Heatmap)\n');
    fprintf('  4. Frequency Analysis (Mag/Phase / Range Profile / Spectrogram)\n');
    fprintf('  5. Bistatic Analysis (Heatmap / Monostatic Cut / Incident Cuts)\n');
    fprintf('  --- Comparison & Export ---\n');
    fprintf('  6. Compare Datasets (Overlay / Difference / Statistics)\n');
    fprintf('  7. Export (CSV / XLSX / PNG / SVG / MAT)\n');
    fprintf('  --- Tools ---\n');
    fprintf('  8. Batch Process (Process all files in a directory)\n');
    fprintf('  0. Exit\n\n');
end

function ok = checkData(data)
    ok = true;
    if isempty(data)
        fprintf(2, '\n  No data loaded. Please load a file first (Option 1).\n');
        ok = false;
    end
end

%% ========================================================================
%  SELECTBISTATICBLOCK  Extract a single incident-angle block for 2D/3D plotting
% ========================================================================
function dataOut = selectBistaticBlock(data, purpose)
    fprintf('\n  --- Select Incident Angle for %s ---\n', purpose);
    fprintf('  Available incident angles: %.1f° ~ %.1f° (%d blocks)\n', ...
            min(data.theta_inc), max(data.theta_inc), data.N_inc);

    % Show the default (center) block
    defaultIdx = data.defaultIncBlock;
    fprintf('  Default: block %d (theta_i = %.1f°)\n', defaultIdx, ...
            data.theta_inc(defaultIdx));

    idxStr = input(sprintf('  Enter block index (1-%d, default=%d): ', ...
                   data.N_inc, defaultIdx), 's');
    idxStr = strtrim(idxStr);
    if isempty(idxStr)
        blkIdx = defaultIdx;
    else
        blkIdx = str2double(idxStr);
    end

    if isnan(blkIdx) || blkIdx < 1 || blkIdx > data.N_inc
        fprintf(2, '  Invalid index. Using default block %d.\n', defaultIdx);
        blkIdx = defaultIdx;
    end

    % Build a non-bistatic data struct compatible with 2D/3D plot functions
    it = data.it;
    ip = data.ip;
    rcsTh = reshape(data.RCS_theta(:, blkIdx), ip, it);
    rcsPh = reshape(data.RCS_phi(:, blkIdx), ip, it);

    dataOut = struct();
    dataOut.type     = 'ffe_single';
    dataOut.filePath = data.filePath;
    dataOut.fileName = sprintf('%s_inci%.0f', data.fileName, data.theta_inc(blkIdx));
    dataOut.theta    = data.theta;
    dataOut.phi      = data.phi;
    dataOut.Sth      = 10 * log10(max(rcsTh, 1e-30));
    dataOut.Sph      = 10 * log10(max(rcsPh, 1e-30));
    dataOut.ip       = ip;
    dataOut.it       = it;
    dataOut.is1D     = (ip == 1) || (it == 1);
    dataOut.is2D     = (ip > 1) && (it > 1);
    dataOut.N_angles = ip * it;
    dataOut.S_complex = [];
    dataOut.freq     = data.freq;
    dataOut.freq_ghz = data.freq_ghz;
    dataOut.N_f      = 1;
    dataOut.hasFreq  = false;
    dataOut.bistatic = false;
    dataOut.param    = sprintf('%s\n  Incident angle block: %d/%d (theta_i = %.1f°)', ...
                               data.param, blkIdx, data.N_inc, data.theta_inc(blkIdx));
    dataOut.paramStruct = data.paramStruct;

    fprintf('  Selected: theta_i = %.1f°\n', data.theta_inc(blkIdx));
end

%% ========================================================================
%  OPTION MENUS
% ========================================================================
function opt = getPlotOptions()
    fprintf('\n  --- 2D Plot Options ---\n');
    comp = input('  Component [theta/phi/both] (default: theta): ', 's');
    comp = strtrim(comp); if isempty(comp), comp = 'theta'; end

    sm = input('  Smoothing window (0=none, default: 0): ', 's');
    sm = strtrim(sm); if isempty(sm), sm = 0; else sm = str2double(sm); end

    bw = input('  Show -3dB beamwidth? [y/n] (default: y): ', 's');
    bw = strtrim(bw); if isempty(bw), bw = true; else bw = strcmpi(bw, 'y'); end

    pk = input('  Show peak markers? [y/n] (default: y): ', 's');
    pk = strtrim(pk); if isempty(pk), pk = true; else pk = strcmpi(pk, 'y'); end

    sv = input('  Save figures? [y/n] (default: n): ', 's');
    sv = strtrim(sv); if isempty(sv), sv = false; else sv = strcmpi(sv, 'y'); end

    opt = struct('component', comp, 'smoothWin', sm, ...
                 'showBeamwidth', bw, 'showPeaks', pk, 'saveFigs', sv);
end

function opt = get3DOptions()
    fprintf('\n  --- 3D Plot Options ---\n');
    comp = input('  Component [theta/phi/both] (default: theta): ', 's');
    comp = strtrim(comp); if isempty(comp), comp = 'theta'; end

    pt = input('  Plot type [all/spherical/mesh/heatmap] (default: all): ', 's');
    pt = strtrim(pt); if isempty(pt), pt = 'all'; end

    cm = input('  Colormap [jet/parula/hot/turbo] (default: jet): ', 's');
    cm = strtrim(cm); if isempty(cm), cm = 'jet'; end

    sv = input('  Save figures? [y/n] (default: n): ', 's');
    sv = strtrim(sv); if isempty(sv), sv = false; else sv = strcmpi(sv, 'y'); end

    opt = struct('component', comp, 'plotType', pt, 'colormap', cm, 'saveFigs', sv);
end

function opt = getFreqOptions(data)
    fprintf('\n  --- Frequency Analysis Options ---\n');
    fprintf('  Total observation angles: %d\n', data.N_angles);
    angIdx = input(sprintf('  Angle index (1-%d, default=center %d): ', ...
                   data.N_angles, round(data.N_angles/2)), 's');
    angIdx = strtrim(angIdx);
    if isempty(angIdx)
        angIdx = round(data.N_angles / 2);
    else
        angIdx = str2double(angIdx);
    end

    wf = input('  Window for IFFT [rect/hamming/hann/blackman] (default: rect): ', 's');
    wf = strtrim(wf); if isempty(wf), wf = 'rect'; end

    zp = input('  Zero-padding factor (default: 4): ', 's');
    zp = strtrim(zp); if isempty(zp), zp = 4; else zp = str2double(zp); end

    sv = input('  Save figures? [y/n] (default: n): ', 's');
    sv = strtrim(sv); if isempty(sv), sv = false; else sv = strcmpi(sv, 'y'); end

    opt = struct('angleIndex', angIdx, 'windowFunc', wf, ...
                 'zeroPad', zp, 'saveFigs', sv);
end

function opt = getBistaticOptions(data)
    fprintf('\n  --- Bistatic Analysis Options ---\n');
    comp = input('  Component [total/theta/phi] (default: total): ', 's');
    comp = strtrim(comp); if isempty(comp), comp = 'total'; end

    fprintf('  Plot types: all, heatmap, monostatic, cuts, forward_back, 3d\n');
    pt = input('  Plot type (default: all): ', 's');
    pt = strtrim(pt); if isempty(pt), pt = 'all'; end

    cm = input('  Colormap [jet/parula/hot/turbo] (default: jet): ', 's');
    cm = strtrim(cm); if isempty(cm), cm = 'jet'; end

    sv = input('  Save figures? [y/n] (default: n): ', 's');
    sv = strtrim(sv); if isempty(sv), sv = false; else sv = strcmpi(sv, 'y'); end

    opt = struct('component', comp, 'plotType', pt, 'colormap', cm, 'saveFigs', sv);
end

function opt = getCompareOptions()
    fprintf('\n  --- Compare Options ---\n');
    comp = input('  Component [theta/phi/both] (default: theta): ', 's');
    comp = strtrim(comp); if isempty(comp), comp = 'theta'; end

    nm = input('  Normalize (align peaks)? [y/n] (default: n): ', 's');
    nm = strtrim(nm); if isempty(nm), nm = false; else nm = strcmpi(nm, 'y'); end

    sv = input('  Save figures? [y/n] (default: n): ', 's');
    sv = strtrim(sv); if isempty(sv), sv = false; else sv = strcmpi(sv, 'y'); end

    opt = struct('component', comp, 'normalize', nm, 'saveFigs', sv);
end

function fmt = getExportFormat()
    fprintf('\n  --- Export Format ---\n');
    fprintf('  Available: csv, xlsx, mat, png, svg, eps, fig, all\n');
    fmt = input('  Format(s) [comma-separated] (default: csv): ', 's');
    fmt = strtrim(fmt);
    if isempty(fmt)
        fmt = 'csv';
    elseif contains(fmt, ',')
        fmt = strsplit(fmt, ',');
        fmt = strtrim(fmt);
    end
end

function opt = getExportOptions()
    svDir = input('  Export directory (default: ./ffe_export/): ', 's');
    svDir = strtrim(svDir);
    if isempty(svDir)
        svDir = './ffe_export';
    end

    exportFigs = input('  Export open figures too? [y/n] (default: n): ', 's');
    exportFigs = strtrim(exportFigs);
    if isempty(exportFigs)
        exportFigs = false;
    else
        exportFigs = strcmpi(exportFigs, 'y');
    end

    opt = struct('saveDir', svDir, 'exportFigs', exportFigs);
end

%% ========================================================================
%  BATCHPROCESS  Process all supported files in a directory
% ========================================================================
function batchProcess(dirPath)
    fprintf('\n========================================\n');
    fprintf('  Batch Processing: %s\n', dirPath);
    fprintf('========================================\n\n');

    % Find all supported files
    ffeFiles = dir(fullfile(dirPath, '*.ffe'));
    datFiles = dir(fullfile(dirPath, '*.dat'));
    matFiles = dir(fullfile(dirPath, '**', '*.mat'));

    allFiles = [ ...
        arrayfun(@(f) fullfile(f.folder, f.name), ffeFiles, 'UniformOutput', false), ...
        arrayfun(@(f) fullfile(f.folder, f.name), datFiles, 'UniformOutput', false), ...
        arrayfun(@(f) fullfile(f.folder, f.name), matFiles, 'UniformOutput', false) ...
    ];

    % Exclude non-result .mat files
    allFiles = allFiles(cellfun(@(f) ...
        contains(f, {'wideband_scattering', 'range_profile', 'sar_imaging', ...
                     'scatter3d', 'sar3d', '.ffe', '.dat'}), allFiles));

    % Also include .ffe files that were filtered out
    for k = 1:length(ffeFiles)
        fpath = fullfile(ffeFiles(k).folder, ffeFiles(k).name);
        if ~any(strcmp(allFiles, fpath))
            allFiles{end+1} = fpath; %#ok<AGROW>
        end
    end

    if isempty(allFiles)
        fprintf('  No result files found in: %s\n', dirPath);
        fprintf('  Looking for: *.ffe, *.dat, wideband_scattering_*.mat, etc.\n');
        return;
    end

    fprintf('  Found %d result files to process.\n\n', length(allFiles));

    % Create export directory
    exportDir = fullfile(dirPath, 'ffe_batch_export');
    if ~exist(exportDir, 'dir')
        mkdir(exportDir);
    end

    for i = 1:length(allFiles)
        filePath = allFiles{i};
        [~, name, ext] = fileparts(filePath);
        fprintf('[%d/%d] %s%s\n', i, length(allFiles), name, ext);

        try
            data = ffe_load_data(filePath);

            if isfield(data, 'bistatic') && data.bistatic
                % Bistatic .ffe: generate bistatic plots
                ffe_plot_bistatic(data, struct('saveFigs', true, ...
                    'saveDir', exportDir, 'figPrefix', name));
                % Also plot center-incident-angle 2D cut
                data2d = selectBistaticBlockAuto(data);
                ffe_plot_rcs_2d(data2d, struct('saveFigs', true, ...
                    'saveDir', exportDir, 'figPrefix', name));
            else
                if data.is1D || data.is2D
                    ffe_plot_rcs_2d(data, struct('saveFigs', true, ...
                        'saveDir', exportDir, 'figPrefix', name));
                end
                if data.is2D
                    ffe_plot_rcs_3d(data, struct('saveFigs', true, ...
                        'saveDir', exportDir, 'figPrefix', name));
                end
                if data.hasFreq
                    ffe_plot_frequency(data, struct('saveFigs', true, ...
                        'saveDir', exportDir, 'figPrefix', name));
                end
            end

            % Export CSV data
            ffe_export(data, 'csv', struct('saveDir', exportDir, ...
                'exportData', true, 'exportFigs', false));

            close all;

        catch ME
            fprintf(2, '  ERROR: %s\n', ME.message);
        end

        fprintf('\n');
    end

    fprintf('========================================\n');
    fprintf('  Batch processing complete.\n');
    fprintf('  Output: %s\n', exportDir);
    fprintf('========================================\n');
end

%% ========================================================================
function dataOut = selectBistaticBlockAuto(data)
    % Auto-select the default (center) incident angle block
    blkIdx = data.defaultIncBlock;
    it = data.it; ip = data.ip;
    rcsTh = reshape(data.RCS_theta(:, blkIdx), ip, it);
    rcsPh = reshape(data.RCS_phi(:, blkIdx), ip, it);

    dataOut = struct();
    dataOut.type     = 'ffe_single';
    dataOut.filePath = data.filePath;
    dataOut.fileName = data.fileName;
    dataOut.theta    = data.theta;
    dataOut.phi      = data.phi;
    dataOut.Sth      = 10 * log10(max(rcsTh, 1e-30));
    dataOut.Sph      = 10 * log10(max(rcsPh, 1e-30));
    dataOut.ip       = ip;
    dataOut.it       = it;
    dataOut.is1D     = (ip == 1) || (it == 1);
    dataOut.is2D     = (ip > 1) && (it > 1);
    dataOut.N_angles = ip * it;
    dataOut.S_complex = [];
    dataOut.freq     = data.freq;
    dataOut.freq_ghz = data.freq_ghz;
    dataOut.N_f      = 1;
    dataOut.hasFreq  = false;
    dataOut.bistatic = false;
    dataOut.param    = data.param;
    dataOut.paramStruct = data.paramStruct;
end
