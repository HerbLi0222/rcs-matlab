function [resultDir, nowStr] = createResultDir(programType)
% CREATERESULTDIR  Create a timestamped result subdirectory under results/
%
%   Creates a directory: results/<programType>_<YYYYmmddHHMMSS>/
%   and returns the full path along with the timestamp string.
%
%   Input:
%       programType - string, program name (e.g., 'bistatic', 'sar_imaging')
%   Output:
%       resultDir   - full path to the created result subdirectory
%       nowStr      - timestamp string (YYYYmmddHHMMSS)
%
%   Usage:
%       [resultDir, nowStr] = createResultDir('bistatic');
%       % resultDir = 'results/bistatic_20260618120000'
%
%   See also: fullfile, datestr

    nowStr = datestr(now, 'yyyymmddHHMMSS');
    dirName = sprintf('%s_%s', programType, nowStr);
    resultDir = fullfile('results', dirName);

    if ~exist(resultDir, 'dir')
        [status, msg] = mkdir(resultDir);
        if ~status
            error('createResultDir:FailedToCreate', ...
                  'Failed to create result directory "%s": %s', resultDir, msg);
        end
    end

    fprintf('  Result directory: %s\n', resultDir);
end
