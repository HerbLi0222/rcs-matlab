function [Et2, Ep2] = biIncidentFieldSphericalCoordinates(cpi2, cti2, sti2, spi2, e2)
    % BIINCIDENTFIELDSPHERICALCOORDINATES Convert incident field to local spherical (bistatic).
    %
    %   Uses pre-computed trig values from the incidence local angles.
    %
    %   Input:
    %       cpi2, cti2, sti2, spi2 - trig functions of local incident angles
    %       e2                    - local Cartesian field vector [3x1]
    %   Output:
    %       Et2 - theta component of incident field
    %       Ep2 - phi component of incident field

    Et2 = e2(1) * cti2 * cpi2 + e2(2) * cti2 * spi2 - e2(3) * sti2;
    Ep2 = -e2(1) * spi2 + e2(2) * cpi2;
end
