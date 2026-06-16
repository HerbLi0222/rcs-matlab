function [perp, para] = reflectionCoefficients(rs, index, th2, thri, phrii, alpha, beta, freq, matrl)
    % REFLECTIONCOEFFICIENTS Compute reflection coefficients for a triangle facet.
    %
    %   Uses either simple resistivity-based model (PEC/finite conductivity)
    %   or detailed material properties from the matrl file.
    %
    %   Input:
    %       rs    - resistivity value for this facet
    %       index - facet index
    %       th2   - local observation theta angle
    %       thri  - incidence theta angle (unused for simple model)
    %       phrii - incidence phi angle (unused for simple model)
    %       alpha - facet azimuth angle
    %       beta  - facet elevation angle
    %       freq  - radar frequency (Hz)
    %       matrl - material properties cell array (may be empty)
    %   Output:
    %       perp - perpendicular (TE) reflection coefficient
    %       para - parallel (TM) reflection coefficient

    MATERIALESPECIFICO = 1;

    if rs == MATERIALESPECIFICO
        % Use detailed material properties
        [perp, para] = getReflCoeffFromMatrl(thri, phrii, alpha, beta, freq, matrl{index});
    else
        % Simple resistivity-based model
        % Local TE polarization
        perp = -1 / (2 * rs * cos(th2) + 1);

        % Local TM polarization
        para = 0;
        if (2 * rs + cos(th2)) ~= 0
            para = -cos(th2) / (2 * rs + cos(th2));
        end
    end
end
