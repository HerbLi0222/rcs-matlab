function [Area, alpha, beta, N, d, ip, it] = calculateValues(pstart, pstop, delp, tstart, tstop, delt, ntria, rad)
    % CALCULATEVALUES Pre-allocate arrays and compute loop bounds for monostatic RCS.
    %
    %   Output:
    %       Area  - (ntria x 1) triangle areas
    %       alpha - (ntria x 1) azimuth angles of normals
    %       beta  - (ntria x 1) elevation angles of normals
    %       N     - (ntria x 3) unit normal vectors
    %       d     - (ntria x 3) edge lengths
    %       ip    - number of phi steps
    %       it    - number of theta steps

    % Calculate number of phi steps
    if delp == 0
        ip = round((pstop - pstart) + 1);
    else
        ip = round((pstop - pstart) / delp + 1);
    end

    % Calculate number of theta steps
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
