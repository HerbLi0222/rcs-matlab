function matrl = getMaterialPropertiesFromFile(filename)
    % GETMATERIALPROPERTIESFROMFILE Parse material properties from a file.
    %
    %   Input:
    %       filename - path to material file or file ID
    %   Output:
    %       matrl - parsed material properties list

    % Open file
    fid = fopen(filename, 'r');
    if fid == -1
        error('getMaterialPropertiesFromFile:FileNotFound', ...
              'Cannot open file: %s', filename);
    end

    % Read all lines
    materialTextList = {};
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line) && ~isempty(strtrim(line))
            materialTextList{end+1} = line; %#ok<AGROW>
        end
    end
    fclose(fid);

    % Convert text to material list
    matrl = convertMaterialTextlistToList(materialTextList);
end
