function coordinatesData = extractCoordinatesData(Rs)
    % EXTRACTCOORDINATESDATA Read and process coordinate/facet data for RCS simulation.
    %
    %   Input:
    %       Rs - resistivity value
    %   Output:
    %       coordinatesData - cell array containing all geometry data:
    %         {x, y, z, xpts, ypts, zpts, nverts, nfc, node1, node2, node3,
    %          iflag, ilum, Rs, ntria, vind, r}

    % Read coordinates
    [x, y, z, xpts, ypts, zpts, nverts] = readCoordinates();

    % Read facets
    [nfc, node1, node2, node3, iflag, ilum, Rs, ntria] = readFacets(Rs);

    % Build vertex index table
    vind = createVind(node1, node2, node3);

    % Build position vectors
    r = calculateR(x, y, z, nverts);

    % Pack into cell array
    coordinatesData = {x, y, z, xpts, ypts, zpts, nverts, nfc, ...
                       node1, node2, node3, iflag, ilum, Rs, ntria, vind, r};
end
