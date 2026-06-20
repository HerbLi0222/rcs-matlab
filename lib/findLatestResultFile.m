function resultFile = findLatestResultFile(pattern)
% FINDLATESTRESULTFILE  Find the most recent result file matching a pattern.
%
%   Recursively searches results/ and all its subdirectories for files
%   matching the given pattern, and returns the full path to the most
%   recently modified file.
%
%   Uses MATLAB's built-in recursive wildcard ('**') on R2016b+.
%   Falls back to flat results/ search on older releases for compatibility.
%
%   Input:
%       pattern    - file pattern string (e.g., 'wideband_scattering_*.mat')
%   Output:
%       resultFile - full path to the newest matching file, or '' if none found
%
%   Usage:
%       f = findLatestResultFile('wideband_scattering_*.mat');
%       if isempty(f)
%           error('No wideband_scattering_*.mat found.');
%       end
%       data = load(f);

    % ---- Try recursive search first (R2016b+) ----
    files = dir(fullfile('results', '**', pattern));

    % ---- Fallback: flat search for older MATLAB ----
    if isempty(files)
        files = dir(fullfile('results', pattern));
    end

    if isempty(files)
        resultFile = '';
        return;
    end

    % ---- Sort by modification date, newest first ----
    [~, idx] = sort([files.datenum], 'descend');
    newest = files(idx(1));

    % ---- Construct full path ----
    % On R2016b+ with '**', the struct has a .folder field.
    % On older MATLAB (flat search), use 'results' as the folder.
    if isfield(newest, 'folder') && ~isempty(newest.folder)
        resultFile = fullfile(newest.folder, newest.name);
    else
        resultFile = fullfile('results', newest.name);
    end
end
