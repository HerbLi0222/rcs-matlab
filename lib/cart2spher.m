function sphericalVector = cart2spher(cartVector)
    % CART2SPHER Convert Cartesian coordinates to spherical coordinates.
    %
    %   Input:
    %       cartVector - [x, y, z] (3-element vector)
    %   Output:
    %       sphericalVector - [R, theta, phi] (3-element vector)

    x = cartVector(1);
    y = cartVector(2);
    z = cartVector(3);

    R = sqrt(x^2 + y^2 + z^2);
    theta = atan2(sqrt(x^2 + y^2), z);
    phi = atan2(y, x);

    sphericalVector = [R; theta; phi];
end
