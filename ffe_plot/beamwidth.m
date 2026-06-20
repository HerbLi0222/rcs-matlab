% BEAMWIDTH  Compute the -3dB (or custom level) beamwidth of an RCS pattern
%   计算RCS方向图的波束宽度
%
%   Input:
%       data     - 1D array of RCS values (dBsm)
%       angle    - 1D array of corresponding angles (deg)
%       level_db - (optional) dB-down level; default = 3 for -3dB
%   Output:
%       bw       - beamwidth in degrees (NaN if not found)
%       leftIdx  - fractional index of left -dB crossing
%       rightIdx - fractional index of right -dB crossing
%
%   Usage:
%     >> [bw, li, ri] = beamwidth(rcs_dB, theta_deg, 3);
%
%   See also: find_peaks_rcs, smooth_rcs

function [bw, leftIdx, rightIdx] = beamwidth(data, angle, level_db)
    if nargin < 3 || isempty(level_db), level_db = 3; end

    data = data(:)';
    angle = angle(:)';

    [peakVal, peakIdx] = max(data);
    threshold = peakVal - level_db;

    % Search left from peak for first crossing below threshold
    leftIdx = NaN;
    for i = peakIdx:-1:2
        if data(i) >= threshold && data(i-1) < threshold
            frac = (threshold - data(i-1)) / (data(i) - data(i-1));
            leftIdx = i - 1 + frac;
            break;
        end
    end

    % Search right from peak
    rightIdx = NaN;
    N = length(data);
    for i = peakIdx:N-1
        if data(i) >= threshold && data(i+1) < threshold
            frac = (threshold - data(i)) / (data(i+1) - data(i));
            rightIdx = i + frac;
            break;
        end
    end

    if isnan(leftIdx) || isnan(rightIdx)
        bw = NaN;
    else
        leftAngle  = interp1(1:N, angle, leftIdx,  'linear', 'extrap');
        rightAngle = interp1(1:N, angle, rightIdx, 'linear', 'extrap');
        bw = abs(rightAngle - leftAngle);
    end
end
