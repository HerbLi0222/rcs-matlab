function [th2, phi2] = sphericalAngles(u2, v2, w2)
    % SPHERICALANGLES Compute spherical angles from local direction cosines.
    %
    %   Input:
    %       u2, v2, w2 - local direction cosine components
    %   Output:
    %       th2  - local theta angle (radians)
    %       phi2 - local phi angle (radians)

    th2 = asin(sqrt(u2^2 + v2^2) * sign(w2));
    phi2 = atan2(v2, u2 + 1e-10);

    if v2 == 0 && u2 + 1e-10 == 0
        phi2 = 0;
    end
end
