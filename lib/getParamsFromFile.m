function paramList = getParamsFromFile(method)
    % GETPARAMSFROMFILE Read simulation parameters from input data file.
    %
    %   Input:
    %       method - 'monostatic' or 'bistatic'
    %   Output:
    %       paramList - cell array of parameters:
    %         [1]  model_name (string)
    %         [2]  frequency (Hz) - converted from GHz
    %         [3]  correlation distance (m)
    %         [4]  standard deviation (m)
    %         [5]  polarization (0=TM-z, 1=TE-z)
    %         [6]  resistivity (0=PEC, 1=material specific)
    %         [7]  start phi (degrees)
    %         [8]  stop phi (degrees)
    %         [9]  delta phi (degrees)
    %         [10] start theta (degrees)
    %         [11] stop theta (degrees)
    %         [12] delta theta (degrees)
    %         [13+] (bistatic only) theta incidence, phi incidence
    %         [end] matrl file path

    inputDataFile = fullfile('input_files', ['input_data_file_' method '.txt']);

    fid = fopen(inputDataFile, 'r');
    if fid == -1
        error('getParamsFromFile:FileNotFound', 'Cannot open: %s', inputDataFile);
    end

    paramList = {};
    lineNum = 0;

    while ~feof(fid)
        line = strtrim(fgetl(fid));
        lineNum = lineNum + 1;

        % Skip empty lines and comments
        if isempty(line) || line(1) == '#'
            continue;
        end

        % Check if line is numeric or string
        [val, status] = str2num(line); %#ok<ST2NM>
        if status && ~isempty(val)
            paramList{end+1} = val; %#ok<AGROW>
        else
            paramList{end+1} = line; %#ok<AGROW>
        end
    end

    fclose(fid);

    % Convert frequency from GHz to Hz (index 2 = FREQUENCY)
    paramList{2} = paramList{2} * 1e9;

    % Convert STL model to coordinates/facets
    modelFile = fullfile('stl_models', paramList{1});
    stlConverter(modelFile);

    % Handle material file path
    % If resistivity is not material-specific, force default matrl.txt
    if paramList{6} ~= 1  % RESISTIVITY ~= MATERIALESPECIFICO
        paramList{end} = 'matrl.txt';
    else
        % Material-specific case: if path is 'configure', fall back to default
        % (CLI version has no GUI for material configuration)
        if strcmp(paramList{end}, 'configure')
            fprintf('  Warning: Material set to ''configure'' but no GUI available.\n');
            fprintf('  Using default matrl.txt (all PEC). Provide a material file for custom materials.\n');
            paramList{end} = 'matrl.txt';
        end
    end

    fprintf('Parameters loaded for %s simulation.\n', method);
    fprintf('  Model: %s\n', paramList{1});
    fprintf('  Frequency: %.2f GHz\n', paramList{2}/1e9);
end
