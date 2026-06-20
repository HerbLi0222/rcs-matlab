% FFE_LOAD_DATA  Unified data loader for FFE post-processing toolbox
%   统一数据加载器 — 支持 POSTFEKO .ffe / .dat / .mat 格式
%
%   Loads RCS simulation results and returns a standardized data struct
%   for downstream processing and visualization.
%
%   Supported formats:
%     .ffe - POSTFEKO Far-Field Export (单站/双站远场数据)
%            Each file may contain multiple blocks with different incident
%            wave directions, forming a bistatic RCS matrix.
%     .dat - open-rcs-matlab monostatic/bistatic output (theta/phi RCS)
%     .mat - Wideband scattering / HRRP / SAR pipeline output
%
%   Input:
%       filePath - (optional) path to data file; if empty, opens file dialog
%   Output:
%       data     - standardized struct (fields vary by format):
%
%          --- Common fields ---
%           .type       - 'ffe', 'dat', or 'mat'
%           .filePath   - full path to the loaded file
%           .fileName   - short file name (without extension)
%           .theta      - observation theta angle array (ip x it) [deg]
%           .phi        - observation phi angle array (ip x it) [deg]
%           .Sth        - RCS theta component (ip x it) [dBsm]
%           .Sph        - RCS phi component (ip x it) [dBsm]
%           .ip         - number of phi observation points
%           .it         - number of theta observation points
%           .is1D       - true if data is a 1D angular cut
%           .is2D       - true if data is a 2D angular grid
%           .N_angles   - total number of observation angles
%           .param      - parameter string (if available)
%           .paramStruct- parsed parameters (struct)
%           .freq       - frequency [Hz] (single value or array)
%           .freq_ghz   - frequency [GHz]
%
%          --- .ffe (POSTFEKO) specific ---
%           .bistatic   - true (multi-incident-angle data)
%           .theta_inc  - incident theta angles [deg] (1 x N_inc)
%           .phi_inc    - incident phi angles [deg]
%           .N_inc      - number of incident angles
%           .RCS_theta  - linear RCS theta (it x N_inc) [m^2]
%           .RCS_phi    - linear RCS phi (it x N_inc) [m^2]
%           .RCS_total  - linear RCS total (it x N_inc) [m^2]
%           .E_theta    - complex E-field theta (it x N_inc)
%           .E_phi      - complex E-field phi (it x N_inc)
%           .source     - model name from POSTFEKO header
%           .ffeDate    - export date string
%           .ffeSolver  - solver name (e.g., 'POSTFEKO')
%
%          --- .mat specific ---
%           .S_complex  - complex scattering field (N_angles x N_f)
%           .hasFreq    - true if multi-frequency data is present
%           .N_f        - number of frequency points
%
%   Usage:
%     >> data = ffe_load_data;
%     >> data = ffe_load_data('ffe_file/result_rocket_po.ffe');
%     >> data = ffe_load_data('results/temp_20260616220000.dat');
%     >> data = ffe_load_data('results/.../wideband_scattering_*.mat');
%
%   See also: ffe_main, ffe_plot_bistatic, ffe_utils

function data = ffe_load_data(filePath)

    %% ---- Select file if not provided ----
    if nargin < 1 || isempty(filePath)
        [fileName, pathName] = uigetfile( ...
            {'*.ffe;*.dat;*.mat', 'All RCS Files (*.ffe, *.dat, *.mat)'; ...
             '*.ffe', 'POSTFEKO Far-Field Export (*.ffe)'; ...
             '*.dat', 'RCS Text Results (*.dat)'; ...
             '*.mat', 'MATLAB Data Files (*.mat)'}, ...
            'Select RCS result file', './');
        if isequal(fileName, 0)
            error('ffe_load_data:NoFileSelected', 'No file selected. Aborting.');
        end
        filePath = fullfile(pathName, fileName);
    end

    % Verify file exists
    if ~exist(filePath, 'file')
        error('ffe_load_data:FileNotFound', 'Cannot find file: %s', filePath);
    end

    [~, shortName, ext] = fileparts(filePath);
    fprintf('========================================\n');
    fprintf('  FFE Data Loader\n');
    fprintf('========================================\n');
    fprintf('  File: %s%s\n', shortName, ext);

    %% ---- Dispatch by extension ----
    switch lower(ext)
        case '.ffe'
            data = loadFfeFile(filePath, shortName);
        case '.dat'
            data = loadDatFile(filePath, shortName);
        case '.mat'
            data = loadMatFile(filePath, shortName);
        otherwise
            error('ffe_load_data:UnknownFormat', ...
                  'Unsupported file format: %s. Use .ffe, .dat or .mat files.', ext);
    end

    %% ---- Derived fields ----
    if ~isfield(data, 'bistatic') || ~data.bistatic
        data.is1D = (data.ip == 1) || (data.it == 1);
        data.is2D = (data.ip > 1) && (data.it > 1);
    else
        data.is1D = (data.ip == 1) || (data.it == 1);
        data.is2D = (data.ip > 1) && (data.it > 1);
    end
    data.hasFreq = isfield(data, 'N_f') && data.N_f > 1;
    data.N_angles = data.ip * data.it;

    %% ---- Print summary ----
    fprintf('\n  Data Summary:\n');
    fprintf('    Type:       %s\n', data.type);
    if isfield(data, 'bistatic') && data.bistatic
        fprintf('    Mode:       Bistatic (multi-incident-angle)\n');
        fprintf('    Obs Grid:   %d phi × %d theta = %d angles\n', ...
                data.ip, data.it, data.N_angles);
        fprintf('    Inc Angles: %d (%.1f° ~ %.1f°)\n', ...
                data.N_inc, min(data.theta_inc), max(data.theta_inc));
    else
        fprintf('    Grid:       %d phi × %d theta = %d angles\n', ...
                data.ip, data.it, data.N_angles);
        if data.is1D
            fprintf('    Data shape: 1D angular cut\n');
        elseif data.is2D
            fprintf('    Data shape: 2D angular grid\n');
        end
    end
    if isfield(data, 'freq_ghz') && ~isempty(data.freq_ghz)
        if length(data.freq_ghz) > 1
            fprintf('    Frequency:  %.2f - %.2f GHz (%d points)\n', ...
                    min(data.freq_ghz), max(data.freq_ghz), length(data.freq_ghz));
        else
            fprintf('    Frequency:  %.3f GHz\n', data.freq_ghz);
        end
    end
    if ~isempty(data.Sth)
        fprintf('    RCS range:  %.1f - %.1f dBsm\n', ...
                min(data.Sth(:)), max(data.Sth(:)));
    end
    fprintf('========================================\n');
end

%% ========================================================================
%  LOADFFEFILE  Parse POSTFEKO .ffe far-field export file
%
%  Format: Multi-block text file. Each block:
%    ##File Type: Far Field                          (file header, once)
%    #Configuration Name: ...
%    #Request Name: FarField_theta_0_359
%    #Frequency: X.XXXXXXXXXE+XX
%    #No. of Theta Samples: N
%    #No. of Phi Samples: M
%    #Incident Wave Direction: (theta, phi)
%    #No. of Header Lines: 1
%    "Theta" "Phi" "Re(Etheta)" ... "RCS(Total)"    (column header)
%    <N*M data rows>
% ========================================================================
function data = loadFfeFile(filePath, shortName)
    %% ---- Read entire file ----
    fid = fopen(filePath, 'r');
    if fid == -1
        error('ffe_load_data:CannotOpen', 'Cannot open file: %s', filePath);
    end
    lines = {};
    while ~feof(fid)
        lines{end+1} = fgetl(fid); %#ok<AGROW>
    end
    fclose(fid);
    N_lines = length(lines);

    %% ---- Parse file-level header (## and ** lines) ----
    source = '?';
    ffeDate = '?';
    solver = '?';
    for i = 1:min(20, N_lines)
        line = strtrim(lines{i});
        if startsWith(line, '##Source:')
            source = strtrim(strrep(line, '##Source:', ''));
        elseif startsWith(line, '##Date:')
            ffeDate = strtrim(strrep(line, '##Date:', ''));
        elseif startsWith(line, '**')
            solver = strtrim(line);
        end
    end

    %% ---- Find all block start positions (#Configuration Name:) ----
    blockStarts = [];
    for i = 1:N_lines
        line = strtrim(lines{i});
        if startsWith(line, '#Configuration Name:')
            blockStarts(end+1) = i; %#ok<AGROW>
        end
    end
    N_blocks = length(blockStarts);

    if N_blocks == 0
        error('ffe_load_data:NoBlocks', ...
              'No #Configuration Name blocks found in .ffe file.');
    end

    fprintf('  Found %d far-field block(s).\n', N_blocks);

    %% ---- Parse first block to get dimensions ----
    firstBlock = parseFfeBlock(lines, blockStarts(1), ...
        findNextBlockStart(blockStarts, 1, N_lines));
    it = firstBlock.N_theta;
    ip = firstBlock.N_phi;
    N_data_rows = it * ip;

    %% ---- Pre-allocate arrays ----
    freq_hz     = zeros(1, N_blocks);
    theta_inc   = zeros(1, N_blocks);  % incident theta
    phi_inc     = zeros(1, N_blocks);  % incident phi
    RCS_theta   = zeros(N_data_rows, N_blocks);  % linear RCS theta
    RCS_phi     = zeros(N_data_rows, N_blocks);  % linear RCS phi
    RCS_total   = zeros(N_data_rows, N_blocks);  % linear RCS total
    E_theta_re  = zeros(N_data_rows, N_blocks);
    E_theta_im  = zeros(N_data_rows, N_blocks);
    E_phi_re    = zeros(N_data_rows, N_blocks);
    E_phi_im    = zeros(N_data_rows, N_blocks);

    % Observation angles (same for all blocks)
    theta_obs = zeros(N_data_rows, 1);
    phi_obs   = zeros(N_data_rows, 1);

    %% ---- Parse all blocks ----
    configNames = cell(1, N_blocks);
    requestNames = cell(1, N_blocks);

    for b = 1:N_blocks
        startLine = blockStarts(b);
        if b < N_blocks
            endLine = blockStarts(b+1) - 1;
        else
            endLine = N_lines;
        end

        block = parseFfeBlock(lines, startLine, endLine);

        % Store block data
        freq_hz(b)   = block.freq;
        theta_inc(b) = block.inc_theta;
        phi_inc(b)   = block.inc_phi;
        configNames{b} = block.configName;
        requestNames{b} = block.requestName;

        % Verify consistent dimensions
        if block.N_theta ~= it || block.N_phi ~= ip
            warning('ffe_load_data:InconsistentDims', ...
                    ['Block %d has different dimensions (%d×%d) than ' ...
                     'block 1 (%d×%d). Data may be inconsistent.'], ...
                    b, block.N_theta, block.N_phi, it, ip);
        end

        % Store observation angles from first block
        if b == 1
            theta_obs = block.theta_obs;
            phi_obs   = block.phi_obs;
        end

        % Extract data columns
        nData = min(N_data_rows, length(block.theta_obs));
        RCS_theta(1:nData, b)  = block.RCS_theta(1:nData);
        RCS_phi(1:nData, b)    = block.RCS_phi(1:nData);
        RCS_total(1:nData, b)  = block.RCS_total(1:nData);
        E_theta_re(1:nData, b) = block.E_theta_re(1:nData);
        E_theta_im(1:nData, b) = block.E_theta_im(1:nData);
        E_phi_re(1:nData, b)   = block.E_phi_re(1:nData);
        E_phi_im(1:nData, b)   = block.E_phi_im(1:nData);
    end

    %% ---- Detect if frequency or incident angle is the swept variable ----
    uniqueFreq    = unique(freq_hz);
    uniqueThInc   = unique(theta_inc);
    uniquePhiInc  = unique(phi_inc);

    if length(uniqueFreq) > 1
        sweptVar = 'frequency';
    else
        sweptVar = 'incident_angle';
    end

    %% ---- Build observation angle matrices (ip x it) ----
    theta_obs_2d = reshape(theta_obs, ip, it);
    phi_obs_2d   = reshape(phi_obs, ip, it);

    %% ---- Build data struct ----
    data.type     = 'ffe';
    data.filePath = filePath;
    data.fileName = shortName;
    data.source   = source;
    data.ffeDate  = ffeDate;
    data.ffeSolver = solver;

    % Observation grid
    data.theta = theta_obs_2d;  % (ip x it)
    data.phi   = phi_obs_2d;    % (ip x it)
    data.ip    = ip;
    data.it    = it;
    data.N_angles = ip * it;

    % Frequency
    if length(uniqueFreq) == 1
        data.freq     = uniqueFreq;
        data.freq_ghz = uniqueFreq / 1e9;
        data.N_f      = 1;
    else
        data.freq     = uniqueFreq(:)';
        data.freq_ghz = data.freq / 1e9;
        data.N_f      = length(uniqueFreq);
    end

    % Bistatic data
    data.bistatic = true;
    data.theta_inc = theta_inc;  % (1 x N_blocks)
    data.phi_inc   = phi_inc;
    data.N_inc     = N_blocks;
    data.sweptVar  = sweptVar;

    % Linear RCS matrices (N_angles x N_inc)
    data.RCS_theta = RCS_theta;
    data.RCS_phi   = RCS_phi;
    data.RCS_total = RCS_total;

    % Complex E-field
    data.E_theta = complex(E_theta_re, E_theta_im);
    data.E_phi   = complex(E_phi_re, E_phi_im);

    % For compatibility with existing 2D/3D plot functions,
    % Sth/Sph store dBsm. When ip==1 (constant phi), Sth/Sph are
    % built from the center incident angle.
    if N_blocks == 1
        % Single block: use directly
        data.Sth = 10 * log10(max(RCS_theta, 1e-30));
        data.Sph = 10 * log10(max(RCS_phi, 1e-30));
    else
        % Multi-block: Sth/Sph are (ip x it) at a chosen incident angle
        % Default: use the block closest to theta_inc = 0° (monostatic-like)
        [~, centerBlock] = min(abs(theta_inc));
        data.defaultIncBlock = centerBlock;
        rcsTh_1block = reshape(RCS_theta(:, centerBlock), ip, it);
        rcsPh_1block = reshape(RCS_phi(:, centerBlock), ip, it);
        data.Sth = 10 * log10(max(rcsTh_1block, 1e-30));
        data.Sph = 10 * log10(max(rcsPh_1block, 1e-30));
    end

    % S_complex is not in .ffe format
    data.S_complex = [];

    % Build param string
    data.param = sprintf(['POSTFEKO Far-Field Export\n', ...
        '  Source: %s\n  Date: %s\n  Solver: %s\n', ...
        '  Blocks: %d\n  Frequency: %.3f GHz\n', ...
        '  Obs Grid: %d phi × %d theta\n', ...
        '  Incident angles: %.1f° ~ %.1f° (step ~%.1f°)\n', ...
        '  Swept variable: %s\n'], ...
        source, ffeDate, solver, N_blocks, ...
        uniqueFreq(1)/1e9, ip, it, ...
        min(theta_inc), max(theta_inc), ...
        median(diff(unique(theta_inc))), sweptVar);
    data.paramStruct = parseParamString(data.param);
    data.paramStruct.configNames = configNames;
    data.paramStruct.requestNames = requestNames;

    fprintf('  Frequency: %.3f GHz\n', uniqueFreq(1)/1e9);
    fprintf('  Observation: %d theta × %d phi\n', it, ip);
    fprintf('  Incident angles: %d blocks (%.1f° ~ %.1f°)\n', ...
            N_blocks, min(theta_inc), max(theta_inc));
    fprintf('  Swept variable: %s\n', sweptVar);
end

%% ========================================================================
%  PARSEFFEBLOCK  Parse a single far-field block from the .ffe file
% ========================================================================
function block = parseFfeBlock(lines, startLine, endLine)
    block = struct();

    % ---- Extract parameter lines (#...) ----
    configName  = '';
    requestName = '';
    freq        = NaN;
    coordSys    = '';
    N_theta     = 0;
    N_phi       = 0;
    resultType  = '';
    inc_theta   = NaN;
    inc_phi     = NaN;
    nHeaderLines = 1;

    paramEndIdx = startLine;
    for i = startLine:min(endLine, startLine + 30)
        line = strtrim(lines{i});
        if startsWith(line, '#Configuration Name:')
            configName = strtrim(strrep(line, '#Configuration Name:', ''));
        elseif startsWith(line, '#Request Name:')
            requestName = strtrim(strrep(line, '#Request Name:', ''));
        elseif startsWith(line, '#Frequency:')
            freq = sscanf(strrep(line, '#Frequency:', ''), '%f');
        elseif startsWith(line, '#Coordinate System:')
            coordSys = strtrim(strrep(line, '#Coordinate System:', ''));
        elseif startsWith(line, '#No. of Theta Samples:')
            N_theta = sscanf(strrep(line, '#No. of Theta Samples:', ''), '%d');
        elseif startsWith(line, '#No. of Phi Samples:')
            N_phi = sscanf(strrep(line, '#No. of Phi Samples:', ''), '%d');
        elseif startsWith(line, '#Result Type:')
            resultType = strtrim(strrep(line, '#Result Type:', ''));
        elseif startsWith(line, '#Incident Wave Direction:')
            dirStr = strrep(line, '#Incident Wave Direction:', '');
            dirStr = strrep(dirStr, '(', '');
            dirStr = strrep(dirStr, ')', '');
            dirVals = sscanf(dirStr, '%f,%f');
            if length(dirVals) >= 2
                inc_theta = dirVals(1);
                inc_phi   = dirVals(2);
            end
        elseif startsWith(line, '#No. of Header Lines:')
            nHeaderLines = sscanf(strrep(line, '#No. of Header Lines:', ''), '%d');
        elseif ~startsWith(line, '#') && ~isempty(line)
            % First non-# line after parameters is the column header
            paramEndIdx = i - 1;
            break;
        end
        paramEndIdx = i;
    end

    % ---- Find data start (after column header + nHeaderLines) ----
    dataStart = paramEndIdx + 1;  % column header line
    if dataStart <= endLine
        % Check if this line looks like a column header
        if contains(lines{dataStart}, '"Theta"')
            dataStart = dataStart + nHeaderLines;
        end
    end

    % ---- Parse data rows ----
    N_expected = N_theta * N_phi;
    theta_obs = zeros(N_expected, 1);
    phi_obs   = zeros(N_expected, 1);
    E_th_re   = zeros(N_expected, 1);
    E_th_im   = zeros(N_expected, 1);
    E_ph_re   = zeros(N_expected, 1);
    E_ph_im   = zeros(N_expected, 1);
    RCS_th    = zeros(N_expected, 1);
    RCS_ph    = zeros(N_expected, 1);
    RCS_tot   = zeros(N_expected, 1);

    row = 0;
    for i = dataStart:endLine
        if i > length(lines), break; end
        line = strtrim(lines{i});
        if isempty(line), continue; end
        % Stop if we hit the next block header
        if startsWith(line, '#') || startsWith(line, '##') || startsWith(line, '**')
            break;
        end

        nums = sscanf(line, '%f');
        if length(nums) >= 9
            row = row + 1;
            if row > N_expected
                % More rows than expected (safety)
                break;
            end
            theta_obs(row) = nums(1);
            phi_obs(row)   = nums(2);
            E_th_re(row)   = nums(3);
            E_th_im(row)   = nums(4);
            E_ph_re(row)   = nums(5);
            E_ph_im(row)   = nums(6);
            RCS_th(row)    = nums(7);
            RCS_ph(row)    = nums(8);
            RCS_tot(row)   = nums(9);
        end
    end

    % Trim to actual rows read
    if row < N_expected
        theta_obs = theta_obs(1:row);
        phi_obs   = phi_obs(1:row);
        E_th_re   = E_th_re(1:row);
        E_th_im   = E_th_im(1:row);
        E_ph_re   = E_ph_re(1:row);
        E_ph_im   = E_ph_im(1:row);
        RCS_th    = RCS_th(1:row);
        RCS_ph    = RCS_ph(1:row);
        RCS_tot   = RCS_tot(1:row);
    end

    % ---- Fill block struct ----
    block.configName  = configName;
    block.requestName = requestName;
    block.freq        = freq;
    block.coordSys    = coordSys;
    block.N_theta     = N_theta;
    block.N_phi       = N_phi;
    block.resultType  = resultType;
    block.inc_theta   = inc_theta;
    block.inc_phi     = inc_phi;
    block.nHeaderLines = nHeaderLines;
    block.theta_obs   = theta_obs;
    block.phi_obs     = phi_obs;
    block.E_theta_re  = E_th_re;
    block.E_theta_im  = E_th_im;
    block.E_phi_re    = E_ph_re;
    block.E_phi_im    = E_ph_im;
    block.RCS_theta   = RCS_th;
    block.RCS_phi     = RCS_ph;
    block.RCS_total   = RCS_tot;
end

%% ========================================================================
function nextStart = findNextBlockStart(blockStarts, currentIdx, N_lines)
    if currentIdx < length(blockStarts)
        nextStart = blockStarts(currentIdx + 1) - 1;
    else
        nextStart = N_lines;
    end
end

%% ========================================================================
%  LOADDATFILE  Parse .dat text file (monostatic/bistatic output format)
% ========================================================================
function data = loadDatFile(filePath, shortName)
    fid = fopen(filePath, 'r');
    if fid == -1
        error('ffe_load_data:CannotOpen', 'Cannot open file: %s', filePath);
    end

    content = {};
    while ~feof(fid)
        content{end+1} = fgetl(fid); %#ok<AGROW>
    end
    fclose(fid);

    % ---- Locate section markers ----
    paramStart    = find(contains(content, 'Simulation Parameters:'), 1);
    thetaStart    = find(contains(content, 'Theta (deg):'), 1);
    rcsThetaStart = find(contains(content, 'RCS Theta (dBsm):'), 1);
    phiStart      = find(contains(content, 'Phi (deg):'), 1);
    rcsPhiStart   = find(contains(content, 'RCS Phi (dBsm):'), 1);

    % ---- Extract parameter string ----
    param = '';
    if ~isempty(paramStart) && ~isempty(thetaStart)
        for k = (paramStart+1):(thetaStart-1)
            line = strtrim(content{k});
            if ~isempty(line)
                param = [param line newline]; %#ok<AGROW>
            end
        end
    end

    % ---- Parse matrix sections ----
    theta = parseMatrixSection(content, thetaStart+1, rcsThetaStart-1);
    Sth   = parseMatrixSection(content, rcsThetaStart+1, phiStart-1);
    phi   = parseMatrixSection(content, phiStart+1, rcsPhiStart-1);
    Sph   = parseMatrixSection(content, rcsPhiStart+1, length(content));

    if isempty(theta) || isempty(Sth)
        if ~isempty(thetaStart)
            theta = parseMatrixSection(content, thetaStart+1, rcsThetaStart-1);
        end
        if isempty(theta)
            warning('ffe_load_data:EmptyData', 'No data found in file.');
        end
    end

    % ---- Build data struct ----
    data.type       = 'dat';
    data.filePath   = filePath;
    data.fileName   = shortName;
    data.theta      = theta;
    data.phi        = phi;
    data.Sth        = Sth;
    data.Sph        = Sph;
    data.ip         = size(theta, 1);
    data.it         = size(theta, 2);
    data.S_complex  = [];
    data.freq       = [];
    data.freq_ghz   = [];
    data.N_f        = 0;
    data.param      = param;
    data.paramStruct = parseParamString(param);
    data.bistatic   = false;
end

%% ========================================================================
%  LOADMATFILE  Load .mat binary file (wideband scattering / HRRP / SAR)
% ========================================================================
function data = loadMatFile(filePath, shortName)
    matVars = load(filePath);

    % ---- Validate required fields ----
    requiredFields = {'S_complex', 'freq_array', 'theta_array', 'phi_array'};
    missingFields = {};
    for i = 1:length(requiredFields)
        if ~isfield(matVars, requiredFields{i})
            missingFields{end+1} = requiredFields{i}; %#ok<AGROW>
        end
    end

    if ~isempty(missingFields)
        warning('ffe_load_data:MissingFields', ...
                'Missing expected fields in .mat file: %s', strjoin(missingFields, ', '));
    end

    % ---- Extract fields with defaults for missing ----
    S_complex  = getFieldSafe(matVars, 'S_complex', []);
    freq_array = getFieldSafe(matVars, 'freq_array', []);
    theta_array = getFieldSafe(matVars, 'theta_array', []);
    phi_array   = getFieldSafe(matVars, 'phi_array', []);
    param       = getFieldSafe(matVars, 'param', '');
    ip          = getFieldSafe(matVars, 'ip', size(theta_array, 1));
    it          = getFieldSafe(matVars, 'it', size(theta_array, 2));
    N_f         = getFieldSafe(matVars, 'N_f', length(freq_array));
    pol         = getFieldSafe(matVars, 'pol', '?');
    inputModel  = getFieldSafe(matVars, 'inputModel', '?');

    % ---- Compute RCS from complex field if available ----
    if ~isempty(S_complex)
        S_mag = abs(S_complex);
        Sth = mean(S_mag, 2);
        Sth = reshape(Sth, size(theta_array));
        Sth_dB = 10 * log10(max(Sth, 1e-30));
        Sph = Sth_dB;
    else
        Sth_dB = [];
        Sph    = [];
    end

    % ---- Build data struct ----
    data.type       = 'mat';
    data.filePath   = filePath;
    data.fileName   = shortName;
    data.theta      = theta_array;
    data.phi        = phi_array;
    data.Sth        = Sth_dB;
    data.Sph        = Sph;
    data.ip         = ip;
    data.it         = it;
    data.S_complex  = S_complex;
    data.freq       = freq_array(:)';
    data.freq_ghz   = freq_array(:)' / 1e9;
    data.N_f        = N_f;
    data.param      = param;
    data.paramStruct = parseParamString(param);
    data.bistatic   = false;

    % ---- Carry over extra metadata fields ----
    extraFields = {'c', 'B', 'delta_r', 'bbox_center', 'pol', 'inputModel', ...
                   'pstart', 'pstop', 'delp', 'tstart', 'tstop', 'delt', ...
                   'ntria', 'freq_start', 'freq_stop', 'wave_array', ...
                   'range_axis', 'M_dB', 'N_angles_total'};
    for i = 1:length(extraFields)
        fname = extraFields{i};
        if isfield(matVars, fname)
            data.(fname) = matVars.(fname);
        end
    end

    fprintf('  .mat file fields: %s\n', strjoin(fieldnames(matVars)', ', '));
end

%% ========================================================================
%  PARSEMATRIXSECTION  Parse a numeric matrix from the .dat file
% ========================================================================
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

%% ========================================================================
%  GETFIELDSAFE  Get a struct field with a default if missing
% ========================================================================
function val = getFieldSafe(s, fieldName, defaultVal)
    if isfield(s, fieldName)
        val = s.(fieldName);
    else
        val = defaultVal;
    end
end

%% ========================================================================
%  PARSEPARAMSTRING  Parse the parameter text block into key-value struct
% ========================================================================
function ps = parseParamString(param)
    ps = struct();
    if isempty(param), return; end

    lines = strsplit(param, newline);
    for i = 1:length(lines)
        line = strtrim(lines{i});
        if isempty(line), continue; end

        colonPos = strfind(line, ':');
        if ~isempty(colonPos)
            key = strtrim(line(1:colonPos(1)-1));
            val = strtrim(line(colonPos(1)+1:end));
            key = matlab.lang.makeValidName(key);
            numVal = str2double(val);
            if ~isnan(numVal)
                ps.(key) = numVal;
            else
                ps.(key) = val;
            end
        end
    end
end
