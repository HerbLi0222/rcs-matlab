function e0 = biIncidentFieldCartesian(uu, vv, ww, cpi, spi, Et, Ep, e0)
    % BIINCIDENTFIELDCARTESIAN Compute incident field in Cartesian (bistatic version).
    %
    %   Uses incident azimuth cos/sin instead of observation azimuth.
    %
    %   Input:
    %       uu, vv, ww - theta-derivative direction cosines
    %       cpi, spi   - cos(phi_i), sin(phi_i) for incidence
    %       Et, Ep     - polarization components
    %       e0         - pre-allocated output array (3x1)
    %   Output:
    %       e0 - incident field vector [Ex; Ey; Ez]

    e0(1) = uu * Et - spi * Ep;
    e0(2) = vv * Et + cpi * Ep;
    e0(3) = ww * Et;
end
