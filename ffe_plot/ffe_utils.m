% FFE_UTILS  Simple utility functions for FFE post-processing toolbox
%
%   Collection of small helper functions. Larger utilities (find_peaks_rcs,
%   beamwidth, smooth_rcs) have been moved to their own files.
%
%   Functions:
%     dB(x)               - Convert linear magnitude to dB scale
%     linear_db(x)        - Convert dB scale to linear magnitude
%     setup_figure        - Create a consistently-styled figure window
%     print_progress      - Display progress bar in command window
%
%   See also: ffe_main, find_peaks_rcs, beamwidth, smooth_rcs

% =========================================================================
%  dB  Convert linear magnitude to decibel scale
% =========================================================================
function y = dB(x)
    floor_val = 1e-30;
    x_safe = max(abs(x), floor_val);
    y = 10 * log10(x_safe);
end

% =========================================================================
%  LINEAR_DB  Convert dB scale back to linear magnitude
% =========================================================================
function y = linear_db(x_db)
    y = 10.^(x_db / 10);
end

% =========================================================================
%  SETUP_FIGURE  Create a figure with consistent styling
% =========================================================================
function fig = setup_figure(name, width, height)
    if nargin < 2 || isempty(width),  width  = 800; end
    if nargin < 3 || isempty(height), height = 600; end

    fig = figure('Name', name, 'NumberTitle', 'off', ...
                 'Position', [100, 100, width, height], ...
                 'Color', 'w');

    set(gca, 'FontSize', 11, 'LineWidth', 1, 'Box', 'on');
    grid on;
    hold on;
end

% =========================================================================
%  PRINT_PROGRESS  Display a text progress bar in the command window
% =========================================================================
function print_progress(i, N, startTime)
    if nargin < 3, startTime = []; end

    pct = 100 * i / N;
    barLen = 40;
    filled = round(barLen * i / N);
    bar = ['[' repmat('=', 1, filled) repmat(' ', 1, barLen - filled) ']'];

    if ~isempty(startTime)
        elapsed = toc(startTime);
        rate = i / max(elapsed, 1e-6);
        eta = (N - i) / max(rate, 1e-6);
        fprintf('\r  %s %5.1f%%  %d/%d  |  elapsed: %.1fs  ETA: %.1fs', ...
                bar, pct, i, N, elapsed, eta);
    else
        fprintf('\r  %s %5.1f%%  %d/%d', bar, pct, i, N);
    end

    if i == N
        fprintf('\n');
    end
end
