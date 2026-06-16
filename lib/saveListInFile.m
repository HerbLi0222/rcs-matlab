function saveListInFile(materialPropertiesList, filePath)
    % SAVELISTINFILE Save material properties list to a text file.
    %
    %   Input:
    %       materialPropertiesList - cell array of material properties
    %       filePath               - output file path

    fid = fopen(filePath, 'w');
    if fid == -1
        error('saveListInFile:CannotWrite', 'Cannot open file for writing: %s', filePath);
    end

    TYPE = 1;
    DESCRIPTION = 2;
    LAYERS_START = 3;

    for r = 1:length(materialPropertiesList)
        row = materialPropertiesList{r};

        % Write type and description
        fprintf(fid, '%s,%s', row{TYPE}, row{DESCRIPTION});

        % Write each layer
        for l = LAYERS_START:length(row)
            layer = row{l};
            for v = 1:length(layer)
                fprintf(fid, ',%g', layer(v));
            end
        end

        fprintf(fid, '\n');
    end

    fclose(fid);
end
