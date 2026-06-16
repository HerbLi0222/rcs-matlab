function [phi, theta, U, V, W, e0, Sth, Sph] = otherVectorComponents(ip, it)
    % OTHERVECTORCOMPONENTS Pre-allocate arrays for angular RCS data.
    %
    %   Output:
    %       phi, theta - (ip x it) angular grid arrays
    %       U, V, W    - (ip x it) direction cosine arrays
    %       e0         - (3 x 1) incident field placeholder
    %       Sth, Sph   - (ip x it) RCS arrays

    phi = zeros(ip, it);
    theta = zeros(ip, it);
    U = zeros(ip, it);
    V = zeros(ip, it);
    W = zeros(ip, it);
    e0 = zeros(3, 1);

    Sth = zeros(ip, it);
    Sph = zeros(ip, it);
end
