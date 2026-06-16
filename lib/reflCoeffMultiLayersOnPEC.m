function [RCperp, RCpar] = reflCoeffMultiLayersOnPEC(thri, phrii, alpha, beta, freq, matrlLine)
    % REFLCOEFFMULTILAYERSONPEC Reflection coefficient for multiple layers on PEC substrate.
    %
    %   Uses the wave matrix (impedance-based) method for N layers backed by PEC.
    %
    %   Input:
    %       thri, phrii - incidence angles (radians)
    %       alpha, beta - facet orientation (radians)
    %       freq        - frequency (Hz)
    %       matrlLine   - material line with N layers
    %   Output:
    %       RCperp, RCpar - reflection coefficients

    LAYERS_START = 3;

    % Rotation to local coordinates
    T21 = rotationTransfMatrix(alpha, beta);
    sphericalVector = spherglobal2local([1; thri; phrii], T21);
    layers = matrlLine(LAYERS_START:end);

    % PEC termination matrix
    PEC = [1, 0; -1, 0];

    WMatrix_par = eye(2);
    WMatrix_perp = eye(2);

    Z0 = 1;
    wave = 3e8 / freq;
    B0 = 2 * pi / wave;
    thinc = sphericalVector(2);

    Z_par = zeros(1, length(layers));
    Z_perp = zeros(1, length(layers));
    gamma_par = zeros(1, length(layers));
    gamma_perp = zeros(1, length(layers));
    tau_par = zeros(1, length(layers));
    tau_perp = zeros(1, length(layers));

    for i = 1:length(layers)
        layer = layers{i};
        erp = layer(1);
        erdp = erp * layer(2);
        erc = erp - 1i * erdp;
        urp = layer(3);
        urdp = layer(4);
        urc = urp - 1i * urdp;
        t = layer(5) * 1e-3;

        Z_par(i) = sqrt(erc/urc - sin(thinc)^2) / (erc/urc * cos(thinc));
        Z_perp(i) = cos(thinc) / sqrt(erc/urc - sin(thinc)^2);

        if i == 1
            gamma_par(i) = (Z_par(i) - Z0) / (Z_par(i) + Z0);
            gamma_perp(i) = (Z_perp(i) - Z0) / (Z_perp(i) + Z0);
        else
            gamma_par(i) = (Z_par(i) - Z_par(i-1)) / (Z_par(i) + Z_par(i-1));
            gamma_perp(i) = (Z_perp(i) - Z_perp(i-1)) / (Z_perp(i) + Z_perp(i-1));
        end

        tau_par(i) = 1 + gamma_par(i);
        tau_perp(i) = 1 + gamma_perp(i);
        phi_calc = B0 * t * sqrt(erc * urc - sin(thinc)^2);

        T_par = [exp(1i * phi_calc), gamma_par(i) * exp(-1i * phi_calc);
                 gamma_par(i) * exp(1i * phi_calc), exp(-1i * phi_calc)];

        WMatrix_par = 1 / tau_par(i) * WMatrix_par * T_par;

        T_perp = [exp(1i * phi_calc), gamma_perp(i) * exp(-1i * phi_calc);
                  gamma_perp(i) * exp(1i * phi_calc), exp(-1i * phi_calc)];

        WMatrix_perp = 1 / tau_perp(i) * WMatrix_perp * T_perp;
    end

    WMatrix_par = WMatrix_par * PEC;
    WMatrix_perp = WMatrix_perp * PEC;

    RCpar = WMatrix_par(2, 1) / WMatrix_par(1, 1);
    RCperp = WMatrix_perp(2, 1) / WMatrix_perp(1, 1);
end
