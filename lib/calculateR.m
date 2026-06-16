function r = calculateR(x, y, z, nverts)
    % CALCULATER Build position vectors for all vertices.
    %
    %   Input:
    %       x, y, z - coordinate vectors
    %       nverts  - number of vertices
    %   Output:
    %       r - (nverts x 3) matrix of position vectors

    r = zeros(nverts, 3);
    for i = 1:nverts
        r(i, :) = [x(i), y(i), z(i)];
    end
end
