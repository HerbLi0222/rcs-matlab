function [u2, v2, w2, T1, T2] = directionCosines(alpha, beta, D0, m)
    % DIRECTIONCOSINES Compute local direction cosines via rotation matrices.
    %
    %   Rotates the global direction vector D0 into local triangle coordinates
    %   using the triangle's orientation angles alpha and beta.
    %
    %   Input:
    %       alpha - azimuth angle array (radians)
    %       beta  - elevation angle array (radians)
    %       D0    - global direction cosine vector [3 x 1]
    %       m     - triangle index
    %
    %   Output:
    %       u2, v2, w2 - local direction cosines
    %       T1, T2     - rotation matrices

    T1 = [cos(alpha(m)),  sin(alpha(m)), 0;
          -sin(alpha(m)), cos(alpha(m)), 0;
          0,              0,             1];

    T2 = [cos(beta(m)), 0, -sin(beta(m));
          0,            1, 0;
          sin(beta(m)), 0, cos(beta(m))];

    D1 = T1 * D0;
    D2 = T2 * D1;

    u2 = D2(1);
    v2 = D2(2);
    w2 = D2(3);
end
