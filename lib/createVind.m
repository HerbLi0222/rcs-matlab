function vind = createVind(node1, node2, node3)
    % CREATEVIND Build vertex index table from node indices.
    %
    %   Input:
    %       node1, node2, node3 - vertex indices for each facet
    %   Output:
    %       vind - (ntria x 3) matrix of vertex indices

    nTriangles = length(node3);
    vind = zeros(nTriangles, 3);

    for i = 1:nTriangles
        vind(i, :) = [node1(i), node2(i), node3(i)];
    end
end
