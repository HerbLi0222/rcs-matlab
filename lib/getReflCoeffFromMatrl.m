function [RCperp, RCpar] = getReflCoeffFromMatrl(thri, phrii, alpha, beta, freq, matrlLine)
    % GETREFLCOEFFFROMMATRL Dispatch reflection coefficient calculation based on material type.
    %
    %   Input:
    %       thri, phrii - incidence spherical angles (radians)
    %       alpha, beta - facet orientation angles (radians)
    %       freq        - radar frequency (Hz)
    %       matrlLine   - cell array: {TYPE, DESCRIPTION, [layer1], ...}
    %   Output:
    %       RCperp - perpendicular polarization reflection coefficient
    %       RCpar  - parallel polarization reflection coefficient

    materialType = matrlLine{1};

    switch materialType
        case 'PEC'
            RCperp = -1;
            RCpar = -1;

        case 'Composito'
            [RCperp, RCpar] = reflCoeffCompo(thri, phrii, alpha, beta, freq, matrlLine);

        case 'Camada de Composito em PEC'
            [RCperp, RCpar] = reflCoeffCompoLayerOnPEC(thri, phrii, alpha, beta, freq, matrlLine);

        case 'Multiplas Camadas'
            [RCperp, RCpar] = reflCoeffMultiLayers(thri, phrii, alpha, beta, freq, matrlLine);

        case 'Multiplas Camadas em PEC'
            [RCperp, RCpar] = reflCoeffMultiLayersOnPEC(thri, phrii, alpha, beta, freq, matrlLine);

        otherwise
            RCperp = 0;
            RCpar = 0;
    end
end
