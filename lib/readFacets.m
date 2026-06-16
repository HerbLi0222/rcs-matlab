function [nfc, node1, node2, node3, iflag, ilum, Rs, ntria] = readFacets(rs)
    % READFACETS Read facet definitions from facets.txt file.
    %
    %   Input:
    %       rs - resistivity value to assign to all facets
    %   Output:
    %       nfc   - facet numbers
    %       node1, node2, node3 - vertex indices for each facet
    %       iflag - illumination flag (always 0)
    %       ilum  - per-facet illumination flag from file
    %       Rs    - resistivity array (all set to rs)
    %       ntria - number of triangles

    fname2 = 'facets.txt';
    facets = dlmread(fname2);

    nfc = facets(:, 1);
    node1 = facets(:, 2);
    node2 = facets(:, 3);
    node3 = facets(:, 4);
    iflag = 0;
    ilum = facets(:, 5);
    Rs = rs * ones(size(facets(:, 5)));
    ntria = length(node3);
end
