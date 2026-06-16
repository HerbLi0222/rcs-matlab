function stlConverter(file_path)
    % STLCONVERTER Convert binary STL file to coordinates.txt and facets.txt
    %   Reads a binary STL file, extracts unique vertices, builds facet
    %   definitions, and writes intermediate text files for RCS computation.
    %
    %   Matches the behavior of stl_module.py from the Python reference.
    %
    %   Input:
    %       file_path - path to .stl file (string)
    %   Output:
    %       (writes coordinates.txt and facets.txt to current directory)

    fid = fopen(file_path, 'rb');
    if fid == -1
        error('stlConverter:FileNotFound', 'Cannot open file: %s', file_path);
    end

    % Skip 80-byte header
    fread(fid, 80, 'uint8=>uint8');

    % Read number of triangles (uint32, little-endian)
    nTriangles = fread(fid, 1, 'uint32=>uint32');

    % Pre-allocate for all vertices (3 per triangle)
    allVertices = zeros(nTriangles * 3, 3);

    for i = 1:nTriangles
        % Read normal vector (3 x float32) - we recompute anyway
        fread(fid, 3, 'float32=>double');

        % Read 3 vertices (9 x float32)
        v1 = fread(fid, 3, 'float32=>double')';
        v2 = fread(fid, 3, 'float32=>double')';
        v3 = fread(fid, 3, 'float32=>double')';

        % Skip 2-byte attribute
        fread(fid, 1, 'uint16=>uint16');

        % Store vertices
        allVertices(3*i-2, :) = v1;
        allVertices(3*i-1, :) = v2;
        allVertices(3*i, :)   = v3;
    end

    fclose(fid);

    % Find unique vertices (matches np.unique with return_index=True)
    [coordinates, ~, ic] = unique(allVertices, 'rows', 'stable');

    nVerts = size(coordinates, 1);

    % Build facet definitions
    facets = zeros(nTriangles, 6);

    for i = 1:nTriangles
        facets(i, 1) = i;          % face number (nfc)
        facets(i, 2) = ic(3*i-2);  % node1
        facets(i, 3) = ic(3*i-1);  % node2
        facets(i, 4) = ic(3*i);    % node3

        % Compute ilum_flag using the Python approach:
        % is_closed_structure = any((coordinates == face[0]).all(axis=1)) and ...
        % Since coordinates contains ALL unique vertices and face vertices
        % came from allVertices which maps to coordinates via ic,
        % the any() check always returns True for valid STL files.
        % So ilum_flag is always 1.
        facets(i, 5) = 1;  % ilum flag

        facets(i, 6) = 0;  % Rs placeholder (set later by readFacets)
    end

    % Write output files
    dlmwrite('coordinates.txt', coordinates, 'delimiter', ' ', 'precision', '%.6f');
    dlmwrite('facets.txt', facets, 'delimiter', ' ', 'precision', '%d');

    fprintf('STL conversion complete: %d vertices, %d facets\n', nVerts, nTriangles);
end
