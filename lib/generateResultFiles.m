function [nowStr, fileName] = generateResultFiles(theta, Sth, phi, Sph, param, ip, resultDir)
    % GENERATERESULTFILES Write RCS simulation results to a .dat file.
    %
    %   Input:
    %       theta     - theta angle array (ip x it)
    %       Sth       - RCS theta-polarization (ip x it)
    %       phi       - phi angle array (ip x it)
    %       Sph       - RCS phi-polarization (ip x it)
    %       param     - parameter string
    %       ip        - number of phi steps
    %       resultDir - (optional) result subdirectory path; defaults to 'results'
    %   Output:
    %       nowStr  - timestamp string
    %       fileName - path to result file

    if nargin < 7 || isempty(resultDir)
        resultDir = 'results';
        if ~exist(resultDir, 'dir'), mkdir(resultDir); end
    end
    nowStr = datestr(now, 'yyyymmddHHMMSS');
    fileName = fullfile(resultDir, ['temp_' nowStr '.dat']);

    fid = fopen(fileName, 'w');

    fprintf(fid, 'RCS SIMULATOR RESULTS %s\n', nowStr);
    fprintf(fid, '\nSimulation Parameters:\n%s', param);
    fprintf(fid, '\nSimulation Results IR Signature:');
    fprintf(fid, '\nTheta (deg):\n');
    for i1 = 1:ip
        fprintf(fid, '%s\n', mat2str(theta(i1, :)));
    end
    fprintf(fid, '\nRCS Theta (dBsm):\n');
    for i1 = 1:ip
        fprintf(fid, '%s\n', mat2str(Sth(i1, :)));
    end
    fprintf(fid, '\nPhi (deg):\n');
    for i1 = 1:ip
        fprintf(fid, '%s\n', mat2str(phi(i1, :)));
    end
    fprintf(fid, '\nRCS Phi (dBsm):\n');
    for i1 = 1:ip
        fprintf(fid, '%s\n', mat2str(Sph(i1, :)));
    end

    fclose(fid);
end
