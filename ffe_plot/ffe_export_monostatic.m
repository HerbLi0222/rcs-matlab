% FFE_EXPORT_MONOSTATIC  Export monostatic RCS data from .ffe files to .mat
%   从.ffe文件提取单站RCS数据并导出为.mat
%
%   Extracts the monostatic RCS diagonal (theta_obs = theta_inc) from
%   POSTFEKO .ffe bistatic data and saves to a .mat file.
%
%   Input:
%       filePaths - string or cell array of .ffe file paths.
%                   If empty, opens file dialog (multi-select).
%       savePath  - (optional) output .mat file path.
%                   Default: './monostatic_export.mat'
%
%   Output .mat file structure:
%       mono(i).name      - source .ffe file name
%       mono(i).freq_GHz  - radar frequency (GHz)
%       mono(i).source    - CAD model name
%       mono(i).method    - solver method (from file name)
%       mono(i).theta_deg - monostatic angle array [deg] (column vector)
%       mono(i).RCS_dBsm  - monostatic RCS [dBsm] (column vector)
%       mono(i).RCS_theta_dBsm - theta-component RCS
%       mono(i).RCS_phi_dBsm   - phi-component RCS
%
%   Usage:
%     >> ffe_export_monostatic('ffe_file/result_rocket_po.ffe');
%     >> ffe_export_monostatic({'*.ffe','*.ffe'}, 'my_mono_data.mat');
%     >> ffe_export_monostatic();  % file dialog
%
%   See also: ffe_load_data, ffe_plot_monostatic

function ffe_export_monostatic(filePaths, savePath)

    %% ---- Handle inputs ----
    if nargin < 1 || isempty(filePaths)
        [fileNames, pathName] = uigetfile( ...
            {'*.ffe', 'POSTFEKO Far-Field Files (*.ffe)'}, ...
            'Select FFE file(s)', 'ffe_file/', 'MultiSelect', 'on');
        if isequal(fileNames, 0)
            disp('Cancelled.'); return;
        end
        if ischar(fileNames), fileNames = {fileNames}; end
        filePaths = cellfun(@(f) fullfile(pathName, f), fileNames, 'UniformOutput', false);
    end

    if ischar(filePaths) || isstring(filePaths)
        filePaths = {char(filePaths)};
    end

    if nargin < 2 || isempty(savePath)
        savePath = 'monostatic_export.mat';
    end

    N = length(filePaths);
    fprintf('=== Monostatic RCS Export ===\n');
    fprintf('  Processing %d file(s)\n', N);

    %% ---- Initialize output struct ----
    mono = struct();
    monoFields = {'name', 'freq_GHz', 'source', 'method', ...
                   'theta_deg', 'RCS_dBsm', 'RCS_theta_dBsm', 'RCS_phi_dBsm'};

    %% ---- Process each file ----
    for i = 1:N
        fprintf('  [%d/%d] %s\n', i, N, filePaths{i});
        data = ffe_load_data(filePaths{i});

        [~, fname] = fileparts(filePaths{i});

        % Detect method from filename
        if contains(fname, '_lpo', 'IgnoreCase', true)
            method = 'LPO';
        elseif contains(fname, '_po', 'IgnoreCase', true)
            method = 'PO';
        elseif contains(fname, '_mom', 'IgnoreCase', true)
            method = 'MoM';
        elseif contains(fname, '_mlfmm', 'IgnoreCase', true)
            method = 'MLFMM';
        else
            method = 'Unknown';
        end

        % Extract monostatic cuts for all 3 components
        [theta, rcsTotal] = extractMono(data.theta, data.theta_inc, data.RCS_total);
        [~,     rcsTheta] = extractMono(data.theta, data.theta_inc, data.RCS_theta);
        [~,     rcsPhi]   = extractMono(data.theta, data.theta_inc, data.RCS_phi);

        mono(i).name           = fname;
        mono(i).freq_GHz       = data.freq_ghz;
        mono(i).source         = data.source;
        mono(i).method         = method;
        mono(i).theta_deg      = theta(:);
        mono(i).RCS_dBsm       = rcsTotal(:);
        mono(i).RCS_theta_dBsm = rcsTheta(:);
        mono(i).RCS_phi_dBsm   = rcsPhi(:);

        fprintf('    %d points | theta: %.1f~%.1f° | RCS peak: %.1f dBsm\n', ...
                length(theta), min(theta), max(theta), max(rcsTotal));
    end

    %% ---- Save ----
    [saveDir, ~] = fileparts(savePath);
    if ~isempty(saveDir) && ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end

    save(savePath, 'mono', '-v7.3');
    fprintf('\n  Saved: %s\n', savePath);
    fprintf('  Variables: mono (1x%d struct)\n', N);
    fprintf('  Fields: %s\n', strjoin(monoFields, ', '));
end

%% ========================================================================
function [angles, rcsDB] = extractMono(theta2d, thetaInc, rcsLinear)
    % Convert to dB first
    rcsDB_full = 10 * log10(max(rcsLinear, 1e-30));

    thetaObs = theta2d(1, :)';
    thInc = thetaInc(:);

    % Find overlapping range
    tMin = max(min(thetaObs), min(thInc));
    tMax = min(max(thetaObs), max(thInc));
    N = min(length(thetaObs), length(thInc));
    angles = linspace(tMin, tMax, N)';

    rcsDB = zeros(N, 1);
    for k = 1:N
        th = angles(k);
        [~, oi] = min(abs(thetaObs - th));
        [~, ii] = min(abs(thInc - th));
        if oi <= size(rcsDB_full, 1) && ii <= size(rcsDB_full, 2)
            rcsDB(k) = rcsDB_full(oi, ii);
        else
            rcsDB(k) = NaN;
        end
    end
    valid = ~isnan(rcsDB);
    angles = angles(valid);
    rcsDB = rcsDB(valid);
end
