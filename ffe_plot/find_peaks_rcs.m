% FIND_PEAKS_RCS  Find local maxima in an RCS pattern
%   查找RCS方向图中的局部峰值
%
%   Input:
%       data      - 1D array of RCS values (dBsm)
%       threshold - (optional) minimum peak prominence in dB; default = 3
%       minSep    - (optional) minimum separation between peaks in samples; default = 3
%   Output:
%       peakVals  - RCS values at peaks (dBsm)
%       peakIdx   - indices of peaks in the data array
%       peakProm  - prominences of each peak (dB)
%
%   Usage:
%     >> [vals, idx, prom] = find_peaks_rcs(rcs_dB, 5, 5);
%
%   See also: beamwidth, smooth_rcs

function [peakVals, peakIdx, peakProm] = find_peaks_rcs(data, threshold, minSep)
    if nargin < 2 || isempty(threshold), threshold = 3; end
    if nargin < 3 || isempty(minSep),    minSep = 3;    end

    data = data(:)';  % ensure row vector
    N = length(data);

    % Find local maxima with minimum separation
    peakIdx = [];
    for i = 2:N-1
        if data(i) > data(i-1) && data(i) >= data(i+1)
            peakIdx(end+1) = i; %#ok<AGROW>
        end
    end

    if isempty(peakIdx)
        peakVals = []; peakProm = [];
        return;
    end

    % Enforce minimum separation (keep the higher peak)
    keep = true(size(peakIdx));
    for i = 1:length(peakIdx)-1
        if keep(i) && (peakIdx(i+1) - peakIdx(i) < minSep)
            if data(peakIdx(i+1)) > data(peakIdx(i))
                keep(i) = false;
            else
                keep(i+1) = false;
            end
        end
    end
    peakIdx = peakIdx(keep);
    peakVals = data(peakIdx);

    % Compute simple prominence (height above nearest lower valley)
    peakProm = zeros(size(peakVals));
    for k = 1:length(peakIdx)
        pi = peakIdx(k);

        % Left valley
        leftVal = peakVals(k);
        for j = pi-1:-1:1
            if data(j) < leftVal, leftVal = data(j); end
            if data(j) > data(pi), break; end
        end
        % Right valley
        rightVal = peakVals(k);
        for j = pi+1:N
            if data(j) < rightVal, rightVal = data(j); end
            if data(j) > data(pi), break; end
        end
        peakProm(k) = peakVals(k) - max(leftVal, rightVal);
    end

    % Filter by prominence threshold
    keepProm = peakProm >= threshold;
    peakIdx  = peakIdx(keepProm);
    peakVals = peakVals(keepProm);
    peakProm = peakProm(keepProm);
end
