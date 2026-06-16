function [sumt, sump, sumdp, sumdt] = calculaCampos(Area, cfac2, corel, th2, wave, Jy2, Ic, ...
    uu, vv, ww, phr, sumt, sump, sumdt, sumdp, m, Jx2, T1, T2)
    % CALCULACAMPOS Compute scattered field contributions from one triangle.
    %
    %   Transforms local scattered fields to global coordinates and
    %   accumulates the total field over all triangles.
    %
    %   Input:
    %       Area  - triangle area array
    %       cfac2 - diffuse scattering factor
    %       corel - correlation length / wavelength
    %       th2   - local theta angle
    %       wave  - wavelength
    %       Jy2, Jx2 - surface current components (local)
    %       Ic    - area integral result
    %       uu, vv, ww - observation direction derivatives
    %       phr   - phi observation angle
    %       sumt, sump, sumdt, sumdp - running sums (updated in place)
    %       m     - triangle index
    %       T1, T2 - rotation matrices
    %
    %   Output:
    %       Updated sumt, sump, sumdt, sumdp

    % Initialize local arrays
    Es0 = zeros(3, 1);
    Es1 = zeros(3, 1);
    Es2 = zeros(3, 1);
    Ed0 = zeros(3, 1);
    Ed1 = zeros(3, 1);
    Ed2 = zeros(3, 1);

    % Add diffuse component
    Edif = cfac2 * Area(m) * (cos(th2)^2) * exp(-(corel * pi * sin(th2) / wave)^2);

    % Scattered field components for triangle m in local coordinates
    Es2(1) = Jx2 * Ic;
    Es2(2) = Jy2 * Ic;
    Es2(3) = 0;

    Ed2(1) = Jx2 * Edif;
    Ed2(2) = Jy2 * Edif;
    Ed2(3) = 0;

    % Transform back to global coordinates
    Es1 = T2.' * Es2;
    Es0 = T1.' * Es1;

    Ed1 = T2.' * Ed2;
    Ed0 = T1.' * Ed1;

    % Project to spherical components
    Ets = uu * Es0(1) + vv * Es0(2) + ww * Es0(3);
    Eps = -sin(phr) * Es0(1) + cos(phr) * Es0(2);

    Etd = uu * Ed0(1) + vv * Ed0(2) + ww * Ed0(3);
    Epd = -sin(phr) * Ed0(1) + cos(phr) * Ed0(2);

    % Accumulate
    sumt = sumt + Ets;
    sumdt = sumdt + abs(Etd);
    sump = sump + Eps;
    sumdp = sumdp + abs(Epd);
end
