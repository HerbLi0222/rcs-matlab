function e0 = incidentFieldCartesian(uu, vv, ww, e0, Et, phr, Ep)
    % INCIDENTFIELDCARTESIAN Compute incident electric field in global Cartesian coordinates.
    %
    %   Input:
    %       uu, vv, ww - theta-derivative direction cosines
    %       e0         - pre-allocated output array (3x1)
    %       Et         - theta polarization component
    %       phr        - phi angle (radians)
    %       Ep         - phi polarization component
    %   Output:
    %       e0 - incident field vector [Ex; Ey; Ez]

    e0(1) = uu * Et - sin(phr) * Ep;
    e0(2) = vv * Et + cos(phr) * Ep;
    e0(3) = ww * Et;
end
