function localVector = spherglobal2local(sphericalVector, T21)
    % SPHERGLOBAL2LOCAL Transform spherical vector from global to local coordinates.
    %
    %   Input:
    %       sphericalVector - [R, theta, phi] in global spherical coordinates
    %       T21             - 3x3 rotation transformation matrix
    %   Output:
    %       localVector     - [R, theta, phi] in local spherical coordinates

    % Convert to Cartesian
    cartVector = spher2cart(sphericalVector);

    % Apply rotation
    cartVector = T21 * cartVector;

    % Convert back to spherical
    localVector = cart2spher(cartVector);
end
