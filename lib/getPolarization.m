function [pol, Et, Ep] = getPolarization(incidentPolarization)
    % GETPOLARIZATION Returns polarization string and incident field components.
    %
    %   Input:
    %       incidentPolarization - 0 for TM-z (theta-polarized), 1 for TE-z (phi-polarized)
    %   Output:
    %       pol - polarization string ('TM-z' or 'TE-z')
    %       Et  - theta component of incident E-field (complex)
    %       Ep  - phi component of incident E-field (complex)

    if incidentPolarization == 0
        % Theta-polarized (TM-z)
        pol = 'TM-z';
        Et = 1 + 1i * 0;
        Ep = 0 + 1i * 0;
    elseif incidentPolarization == 1
        % Phi-polarized (TE-z)
        pol = 'TE-z';
        Et = 0 + 1i * 0;
        Ep = 1 + 1i * 0;
    else
        error('getPolarization:InvalidInput', 'Invalid polarization value: %d', incidentPolarization);
    end
end
