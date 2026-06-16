function [gammapar, gammaperp, thetat, TIR] = reflCoeff(er1, mr1, er2, mr2, thetai)
    % REFLCOEFF Compute Fresnel reflection coefficients at an interface.
    %
    %   Input:
    %       er1, mr1 - relative permittivity and permeability of medium 1
    %       er2, mr2 - relative permittivity and permeability of medium 2
    %       thetai   - incidence angle (radians)
    %   Output:
    %       gammapar  - parallel polarization reflection coefficient
    %       gammaperp - perpendicular polarization reflection coefficient
    %       thetat    - transmission angle (radians)
    %       TIR       - total internal reflection flag (0 or 1)

    m0 = 4 * pi * 1e-7;   % vacuum permeability
    e0 = 8.854e-12;        % vacuum permittivity

    TIR = 0;

    % Snell's law for transmission angle
    sinthetat = sin(thetai) * sqrt(real(er1) * real(mr1) / (real(er2) * real(mr2)));

    % Check for total internal reflection
    if sinthetat > 1
        TIR = 1;
        thetat = pi / 2;  % critical angle
    else
        thetat = asin(sinthetat);
    end

    % Calculate refractive indices
    n1 = sqrt(mr1 * m0 / (er1 * e0));
    n2 = sqrt(mr2 * m0 / (er2 * e0));

    % Fresnel reflection coefficients
    gammaperp = (n2 * cos(thetai) - n1 * cos(thetat)) / (n2 * cos(thetai) + n1 * cos(thetat));
    gammapar  = (n2 * cos(thetat) - n1 * cos(thetai)) / (n2 * cos(thetat) + n1 * cos(thetai));
end
