function T21 = rotationTransfMatrix(alpha, beta)
    % ROTATIONTRANSFMATRIX Compute combined rotation transformation matrix.
    %   T21 = T2 * T1 where T1 rotates about z-axis by alpha and
    %   T2 rotates about y-axis by beta.
    %
    %   Input:
    %       alpha - rotation angle about z-axis (radians)
    %       beta  - rotation angle about y-axis (radians)
    %   Output:
    %       T21 - 3x3 rotation matrix

    T1 = [cos(alpha),  sin(alpha), 0;
          -sin(alpha), cos(alpha), 0;
          0,           0,          1];

    T2 = [cos(beta), 0, -sin(beta);
          0,         1, 0;
          sin(beta), 0, cos(beta)];

    T21 = T2 * T1;
end
