function [x, y, z, xpts, ypts, zpts, nverts] = readCoordinates()
    % READCOORDINATES Read vertex coordinates from coordinates.txt file.
    %
    %   Output:
    %       x, y, z    - separate coordinate vectors (nverts x 1)
    %       xpts, ypts, zpts - copy of x, y, z for compatibility
    %       nverts     - number of vertices

    fname = 'coordinates.txt';
    coordinates = dlmread(fname);

    xpts = coordinates(:, 1);
    ypts = coordinates(:, 2);
    zpts = coordinates(:, 3);

    x = xpts;
    y = ypts;
    z = zpts;

    nverts = length(xpts);
end
