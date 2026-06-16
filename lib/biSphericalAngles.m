function [thi2, phii2, cpi2, spi2, sti2, cti2] = biSphericalAngles(ui2, vi2, wi2)
    % BISPHERICALANGLES Compute spherical angles from local incidence direction cosines.
    %
    %   Bistatic version: returns additional trig function values.
    %
    %   Input:
    %       ui2, vi2, wi2 - local incidence direction cosine components
    %   Output:
    %       thi2  - local theta angle (radians)
    %       phii2 - local phi angle (radians)
    %       cpi2, spi2, sti2, cti2 - trig functions of thi2, phii2

    sti2 = sqrt(ui2^2 + vi2^2) * sign(wi2);
    cti2 = sqrt(1 - sti2^2);
    thi2 = acos(cti2);
    phii2 = atan2(vi2, ui2 + 1e-10);

    if vi2 == 0 && ui2 + 1e-10 == 0
        phii2 = 0;
    end

    cpi2 = cos(phii2);
    spi2 = sin(phii2);
end
