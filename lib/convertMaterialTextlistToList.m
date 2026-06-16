function matrl = convertMaterialTextlistToList(textList)
    % CONVERTMATERIALTEXTLISTTOLIST Convert comma-separated material text to structured list.
    %
    %   Input:
    %       textList - cell array of strings, each line is comma-separated:
    %                  TYPE,DESCRIPTION,layer1_val1,layer1_val2,...,layerN_val5
    %   Output:
    %       matrl - cell array where each row is:
    %               {TYPE, DESCRIPTION, [layer1], [layer2], ...}
    %               Each layer is [erp, lossTan*erp, urp, urdp, thickness_mm]

    matrl = {};
    TYPE = 1;
    DESCRIPTION = 2;
    LAYERS_START = 3;  % 1-based index where layer data starts

    for r = 1:length(textList)
        % Split by comma
        entries = strsplit(strtrim(textList{r}), ',');

        % First two entries are type and description
        formattedEntry = cell(1, 2);
        formattedEntry{TYPE} = entries{TYPE};
        formattedEntry{DESCRIPTION} = entries{DESCRIPTION};

        % Remaining entries are grouped into layers of 5 values each
        layerData = entries(LAYERS_START:end);
        layerIdx = 1;
        currentLayer = [];

        for k = 1:length(layerData)
            currentLayer(end+1) = str2double(layerData{k}); %#ok<AGROW>

            if mod(k, 5) == 0
                formattedEntry{end+1} = currentLayer; %#ok<AGROW>
                currentLayer = [];
                layerIdx = layerIdx + 1;
            end
        end

        % If there's a partial layer (shouldn't happen with valid data)
        if ~isempty(currentLayer)
            formattedEntry{end+1} = currentLayer; %#ok<AGROW>
        end

        matrl{end+1} = formattedEntry; %#ok<AGROW>
    end
end
