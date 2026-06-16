function [N, d, Area, beta, alpha] = productVector(ntria, N, r, d, Area, alpha, beta, vind)
    % PRODUCTVECTOR Compute triangle normals, areas, edge lengths, and orientation angles.
    %
    %   For each triangle, computes:
    %   - Edge vectors A, B, C
    %   - Outward-facing unit normal N
    %   - Edge lengths d
    %   - Triangle area via Heron's formula
    %   - Elevation angle beta (from z-axis)
    %   - Azimuth angle alpha (from x-axis)
    %
    %   Input:
    %       ntria - number of triangles
    %       N     - pre-allocated normal array
    %       r     - vertex position vectors
    %       d     - pre-allocated edge length array
    %       Area  - pre-allocated area array
    %       alpha - pre-allocated azimuth array
    %       beta  - pre-allocated elevation array
    %       vind  - vertex index table
    %
    %   Output:
    %       Updated arrays with computed values.

    for i = 1:ntria
        % Edge vectors
        A = r(vind(i, 2), :) - r(vind(i, 1), :);
        B = r(vind(i, 3), :) - r(vind(i, 2), :);
        C = r(vind(i, 1), :) - r(vind(i, 3), :);

        % Normal (negative cross product for outward-pointing)
        N(i, :) = -cross(B, A);

        % Edge lengths
        d(i, 1) = norm(A);
        d(i, 2) = norm(B);
        d(i, 3) = norm(C);

        % Area via Heron's formula
        ss = 0.5 * sum(d(i, :));
        Area(i) = sqrt(ss * (ss - d(i, 1)) * (ss - d(i, 2)) * (ss - d(i, 3)));

        % Unit normal
        Nn = norm(N(i, :));
        if Nn ~= 0
            N(i, :) = N(i, :) / Nn;
        end

        % Orientation angles
        % beta: angle from z-axis, 0 < beta < 180
        beta(i) = acos(N(i, 3));

        % alpha: azimuth angle from x-axis
        alpha(i) = atan2(N(i, 2), N(i, 1));

        % Convention for (0,0) normal
        if N(i, 1) == 0 && N(i, 2) == 0
            alpha(i) = 0;
        end
    end
end
