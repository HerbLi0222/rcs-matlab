function matrl = getEntrysFromMatrlFile(ntria, matrlpath)
    % GETENTRYSFROMMATRLFILE Read and validate material properties file.
    %
    %   If the file doesn't exist or has wrong number of entries,
    %   auto-generates a default PEC material list for all facets.
    %
    %   Input:
    %       ntria     - expected number of facets
    %       matrlpath - path to material file
    %   Output:
    %       matrl - cell array of material properties per facet

    % Try to read the material file
    matrl = {};
    try
        matrl = getMaterialPropertiesFromFile(matrlpath);
    catch
        fprintf('  Material file ''%s'' not found. Generating default PEC material.\n', matrlpath);
    end

    % Validate number of entries
    if length(matrl) ~= ntria
        if ~isempty(matrl)
            fprintf('  Material file has %d entries but model has %d facets.\n', ...
                    length(matrl), ntria);
        end
        fprintf('  Auto-generating PEC material for all %d facets.\n', ntria);

        % Generate default PEC material for all facets
        matrl = cell(1, ntria);
        for i = 1:ntria
            matrl{i} = {'PEC', 'facet description', [0, 0, 0, 0, 0]};
        end

        % Save to file for future use
        saveListInFile(matrl, matrlpath);
        fprintf('  Saved default material to ''%s''.\n', matrlpath);
    end
end
