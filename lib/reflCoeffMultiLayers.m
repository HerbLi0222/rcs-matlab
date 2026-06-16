function [RCperp, RCpar] = reflCoeffMultiLayers(thri, phrii, alpha, beta, freq, matrlLine)
    % REFLCOEFFMULTILAYERS Reflection coefficient for multiple dielectric layers (air both sides).
    %
    %   Uses cascaded 2x2 transfer matrices for N layers between two air half-spaces.
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
    THETA = 2;
    layers = matrlLine(LAYERS_START:end);

    Mpar = eye(2);
    Mperp = eye(2);

    er = cell(1, length(layers));
    mr = cell(1, length(layers));
    t = zeros(1, length(layers));
    thetat = zeros(1, length(layers));

    for i = 1:length(layers)
        layer = layers{i};
        er{i} = layer(1) - 1i * layer(2) * layer(1);
        mr{i} = layer(3) - 1i * layer(4);
        t(i) = layer(5) * 0.001;

        if i == 1
            [Gpar, Gperp, thetatI, TIR] = reflCoeff(1, 1, er{i}, mr{i}, sphericalVector(THETA));
        else
            [Gpar, Gperp, thetatI, TIR] = reflCoeff(er{i-1}, mr{i-1}, er{i}, mr{i}, thetat(i-1));
        end

        thetat(i) = thetatI;
        v = 3e8 / sqrt(real(er{i}) * real(mr{i}));
        wave = v / freq;
        b1 = 2 * pi / wave;
        phase = b1 * t(i);

        Mpar = Mpar * [exp(1i * phase), Gpar * exp(-1i * phase);
                        Gpar * exp(1i * phase), exp(-1i * phase)];
        Mperp = Mperp * [exp(1i * phase), Gperp * exp(-1i * phase);
                          Gperp * exp(1i * phase), exp(-1i * phase)];
    end

    % Final interface: last layer → air
    [Gpar, Gperp, thetatdum, TIR] = reflCoeff(er{end}, mr{end}, 1, 1, thetat(end));

    Mpar = Mpar * [exp(1i * phase), Gpar * exp(-1i * phase);
                    Gpar * exp(1i * phase), exp(-1i * phase)];
    Mperp = Mperp * [exp(1i * phase), Gperp * exp(-1i * phase);
                      Gperp * exp(1i * phase), exp(-1i * phase)];

    RCpar = Mpar(2, 1) / Mpar(1, 1);
    RCperp = Mperp(2, 1) / Mperp(1, 1);
end
