function [Area, alpha, beta, N, d, ip, it, cpi, spi, sti, cti, ui, vi, wi, D0i, uui, vvi, wwi, Ri] = ...
    biCalculateValues(pstart, pstop, delp, tstart, tstop, delt, ntria, rad, fii, thetai)
    % BICALCULATEVALUES Pre-compute arrays and incidence vectors for bistatic RCS.
    %
    %   Input:
    %       fii    - incident phi angle (degrees)
    %       thetai - incident theta angle (degrees)
    %   Output (in addition to calculateValues outputs):
    %       cpi, spi, sti, cti - trig functions of incidence angles
    %       ui, vi, wi         - incidence direction vector components
    %       D0i                - incidence direction cosine vector
    %       uui, vvi, wwi      - rotated incidence components
    %       Ri                 - incidence radial unit vector

    % Trig functions of incidence angles
    cpi = cos(fii * pi / 180.0);
    spi = sin(fii * pi / 180.0);
    sti = sin(thetai * pi / 180.0);
    cti = cos(thetai * pi / 180.0);

    % Incidence direction vector
    ui = sti * cpi;
    vi = sti * spi;
    wi = cti;
    D0i = [ui; vi; wi];

    % Rotated incidence components
    uui = cti * cpi;
    vvi = cti * spi;
    wwi = -sti;
    Ri = [ui; vi; wi];

    % Number of phi and theta steps (same as monostatic)
    if delp == 0
        ip = round((pstop - pstart) + 1);
    else
        ip = round((pstop - pstart) / delp + 1);
    end

    if delt == 0
        it = round((tstop - tstart) + 1);
    else
        it = round((tstop - tstart) / delt + 1);
    end

    % Pre-allocate arrays
    Area = zeros(ntria, 1);
    alpha = zeros(ntria, 1);
    beta = zeros(ntria, 1);
    N = zeros(ntria, 3);
    d = zeros(ntria, 3);
end
