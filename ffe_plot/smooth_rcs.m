% SMOOTH_RCS  Apply moving-average smoothing to an RCS pattern
%   对RCS数据应用滑动平均平滑
%
%   Input:
%       data       - 1D array of RCS values (dBsm)
%       windowSize - (optional) smoothing window width; default = 5
%   Output:
%       smoothed - smoothed data (same length as input)
%
%   Usage:
%     >> rcs_smooth = smooth_rcs(rcs_dB, 3);
%
%   See also: find_peaks_rcs, beamwidth

function smoothed = smooth_rcs(data, windowSize)
    if nargin < 2 || isempty(windowSize), windowSize = 5; end

    data = data(:);
    N = length(data);

    if windowSize >= N
        smoothed = repmat(mean(data), N, 1);
        return;
    end

    % Use MATLAB's smoothdata if available (R2017a+); otherwise manual
    try
        smoothed = smoothdata(data, 'movmean', windowSize);
    catch
        halfW = floor(windowSize / 2);
        smoothed = zeros(N, 1);
        for i = 1:N
            lo = max(1, i - halfW);
            hi = min(N, i + halfW);
            smoothed(i) = mean(data(lo:hi));
        end
    end
end
