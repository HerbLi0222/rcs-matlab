function [Et2, Ep2] = incidentFieldSphericalCoordinates(th2, e2, phi2)
    % INCIDENTFIELDSPHERICALCOORDINATES Convert incident field to local spherical coordinates.
    %
    %   Input:
    %       th2  - local theta angle (radians)
    %       e2   - local Cartesian field vector [3x1]
    %       phi2 - local phi angle (radians)
    %   Output:
    %       Et2 - theta component of incident field
    %       Ep2 - phi component of incident field

    Et2 = e2(1) * cos(th2) * cos(phi2) + e2(2) * cos(th2) * sin(phi2) - e2(3) * sin(th2);
    Ep2 = -e2(1) * sin(phi2) + e2(2) * cos(phi2);
end
