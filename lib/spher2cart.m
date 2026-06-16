function cartVector = spher2cart(sphericalVector)
    % SPHER2CART Convert spherical coordinates to Cartesian coordinates.
    %
    %   Input:
    %       sphericalVector - [R, theta, phi] (3-element vector)
    %   Output:
    %       cartVector - [x, y, z] (3-element vector)

    R = sphericalVector(1);
    theta = sphericalVector(2);
    phi = sphericalVector(3);

    x = R * sin(theta) * cos(phi);
    y = R * sin(theta) * sin(phi);
    z = R * cos(theta);

    cartVector = [x; y; z];
end
