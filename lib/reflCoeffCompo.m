function [RCperp, RCpar] = reflCoeffCompo(thri, phrii, alpha, beta, freq, matrlLine)
    % REFLCOEFFCOMPO Reflection coefficient for a single composite layer.
    %
    %   Models a single dielectric/magnetic layer on an infinite substrate.
    %   Uses 2x2 transfer matrix method.
    %
    %   Input:
    %       thri, phrii - incidence angles (radians)
    %       alpha, beta - facet orientation (radians)
    %       freq        - frequency (Hz)
    %       matrlLine   - material line with one layer
    %   Output:
    %       RCperp, RCpar - reflection coefficients

    LAYERS_START = 3;
    layer = matrlLine{LAYERS_START};

    er = layer(1) - 1i * layer(2) * layer(1);
    mr = layer(3) - 1i * layer(4);
    t = layer(5) * 0.001;  % mm to m

    % Rotation to local coordinates
    T21 = rotationTransfMatrix(alpha, beta);
    sphericalVector = spherglobal2local([1; thri; phrii], T21);
    THETA = 2;  % index of theta in spherical vector

    % Fresnel coefficients at interface 1 (air → material)
    [G1par, G1perp, thetat, TIR] = reflCoeff(1, 1, er, mr, sphericalVector(THETA));

    % Fresnel coefficients at interface 2 (material → air, sign reversed)
    G2par = -G1par;
    G2perp = -G1perp;

    % Phase shift through layer
    v = 3e8 / sqrt(real(er) * real(mr));
    wave = v / freq;
    b1 = 2 * pi / wave;
    phase = b1 * t;

    % Transfer matrices
    M1par = [exp(1i * phase), G1par * exp(-1i * phase);
             G1par * exp(1i * phase), exp(-1i * phase)];

    M1perp = [exp(1i * phase), G1perp * exp(-1i * phase);
              G1perp * exp(1i * phase), exp(-1i * phase)];

    M2par = [1, G2par; G2par, 1];
    M2perp = [1, G2perp; G2perp, 1];

    Mpar = M1par * M2par;
    Mperp = M1perp * M2perp;

    RCpar = Mpar(2, 1) / Mpar(1, 1);
    RCperp = Mperp(2, 1) / Mperp(1, 1);
end
