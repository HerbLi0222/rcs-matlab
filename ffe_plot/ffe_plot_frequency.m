% FFE_PLOT_FREQUENCY  Frequency-domain analysis for wideband scattering data
%   宽带散射场频域分析与快速距离像查看
%
%   Analyzes complex scattering field data (S_complex) from wideband
%   simulations. Creates:
%     - Magnitude & phase vs frequency at a selected angle
%     - Quick range profile via IFFT (one-angle HRRP)
%     - Frequency-angle spectrogram (|S| vs frequency vs angle)
%     - Frequency-averaged RCS pattern vs angle
%
%   Requires data loaded from a wideband_scattering_*.mat file.
%
%   Input:
%       data    - data struct from ffe_load_data (must have .S_complex and .freq)
%       options - (optional) struct with fields:
%           .angleIndex  - observation angle index for 1D cuts (default = center)
%           .windowFunc  - window for IFFT: 'rect' (default), 'hamming', 'hann',
%                          'blackman', 'kaiser'
%           .zeroPad     - zero-padding factor for IFFT (default = 4)
%           .saveFigs    - true/false to save figures (default = false)
%           .saveDir     - directory for saved figures
%           .figPrefix   - prefix for saved figure filenames
%
%   Usage:
%     >> data = ffe_load_data('results/.../wideband_scattering_*.mat');
%     >> ffe_plot_frequency(data);
%     >> ffe_plot_frequency(data, struct('windowFunc', 'hamming', 'zeroPad', 8));
%
%   See also: ffe_main, ffe_load_data, main_range_profile

function ffe_plot_frequency(data, options)

    %% ---- Validate ----
    if ~data.hasFreq
        error('ffe_plot_frequency:NoFreqData', ...
              ['No frequency-domain data available. ', ...
               'Load a wideband_scattering_*.mat file.']);
    end

    %% ---- Default options ----
    if nargin < 2, options = struct(); end
    opt = setDefaults(options, data);

    %% ---- Convenience variables ----
    S     = data.S_complex;   % (N_angles x N_f)
    freq  = data.freq;        % [Hz]
    fGHz  = data.freq_ghz;    % [GHz]
    theta = data.theta;
    phi   = data.phi;
    N_angles = data.N_angles;
    N_f  = data.N_f;

    c = 3e8;  % speed of light

    %% ---- Select observation angle ----
    angleIdx = opt.angleIndex;
    if angleIdx > N_angles
        warning('ffe_plot_frequency:BadIndex', ...
                'angleIndex %d > N_angles %d. Using center.', angleIdx, N_angles);
        angleIdx = round(N_angles / 2);
    end

    % Get the theta/phi values for this index
    [phiIdx, thetaIdx] = ind2sub([data.ip, data.it], angleIdx);
    thVal = theta(phiIdx, thetaIdx);
    phVal = phi(phiIdx, thetaIdx);

    fprintf('  Frequency Analysis at angle #%d: theta=%.1f°, phi=%.1f°\n', ...
            angleIdx, thVal, phVal);

    %% ---- Figure 1: Magnitude & Phase vs Frequency ----
    plotMagPhase(freq, fGHz, S, angleIdx, data, opt);

    %% ---- Figure 2: Quick Range Profile (IFFT) ----
    plotQuickRangeProfile(freq, S, angleIdx, thVal, phVal, N_f, c, data, opt);

    %% ---- Figure 3: Spectrogram (|S| vs frequency vs angle) ----
    plotSpectrogram(freq, fGHz, S, theta, phi, data, opt);

    %% ---- Figure 4: Frequency-Averaged RCS vs Angle ----
    plotFreqAveragedRCS(S, theta, phi, data, opt);
end

%% ========================================================================
function opt = setDefaults(options, data)
    opt.angleIndex = getOpt(options, 'angleIndex', round(data.N_angles / 2));
    opt.windowFunc = getOpt(options, 'windowFunc', 'rect');
    opt.zeroPad    = getOpt(options, 'zeroPad', 4);
    opt.saveFigs   = getOpt(options, 'saveFigs', false);
    opt.saveDir    = getOpt(options, 'saveDir', fullfile('..', 'results'));
    opt.figPrefix  = getOpt(options, 'figPrefix', data.fileName);
end

function val = getOpt(s, field, defaultVal)
    if isfield(s, field) && ~isempty(s.(field))
        val = s.(field);
    else
        val = defaultVal;
    end
end

%% ========================================================================
%  Figure 1: Magnitude & Phase vs Frequency
% ========================================================================
function plotMagPhase(freq, fGHz, S, angleIdx, data, opt)
    fig = figure('Name', sprintf('FFE - Mag/Phase [%s]', data.fileName), ...
                 'NumberTitle', 'off', ...
                 'Position', [50, 80, 950, 500], 'Color', 'w');

    S_sel = S(angleIdx, :);
    mag_dB = 20 * log10(max(abs(S_sel), 1e-30));
    phase_deg = unwrap(angle(S_sel)) * 180 / pi;

    % ---- Magnitude subplot ----
    subplot(2, 1, 1);
    plot(fGHz, mag_dB, 'b-', 'LineWidth', 1.8);
    xlabel('Frequency (GHz)', 'FontSize', 12);
    ylabel('|S| (dB)', 'FontSize', 12);
    title(sprintf('Scattering Magnitude vs Frequency  (angle #%d)', angleIdx), ...
          'FontSize', 11);
    grid on; box on;
    xlim([min(fGHz) max(fGHz)]);

    % Mark peak response frequency
    [peakMag, peakI] = max(mag_dB);
    hold on;
    plot(fGHz(peakI), peakMag, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    text(fGHz(peakI), peakMag, ...
         sprintf('  %.1f dB @ %.2f GHz', peakMag, fGHz(peakI)), ...
         'FontSize', 9, 'Color', 'r');

    % ---- Phase subplot ----
    subplot(2, 1, 2);
    plot(fGHz, phase_deg, 'r-', 'LineWidth', 1.5);
    xlabel('Frequency (GHz)', 'FontSize', 12);
    ylabel('Phase (deg)', 'FontSize', 12);
    title('Unwrapped Phase vs Frequency', 'FontSize', 11);
    grid on; box on;
    xlim([min(fGHz) max(fGHz)]);

    sgtitle(sprintf('Frequency Response  —  %s', data.fileName), 'FontSize', 12);

    if opt.saveFigs
        saveFig(fig, fullfile(opt.saveDir, [opt.figPrefix '_freq_magphase']), opt);
    end
end

%% ========================================================================
%  Figure 2: Quick Range Profile via IFFT
% ========================================================================
function plotQuickRangeProfile(freq, S, angleIdx, thVal, phVal, N_f, c, data, opt)
    fig = figure('Name', sprintf('FFE - Quick Range Profile [%s]', data.fileName), ...
                 'NumberTitle', 'off', ...
                 'Position', [50, 60, 750, 500], 'Color', 'w');

    S_sel = S(angleIdx, :);

    % ---- Apply window ----
    win = createWindow(N_f, opt.windowFunc);
    S_windowed = S_sel .* win(:)';

    % ---- Zero-padding ----
    N_fft = N_f * opt.zeroPad;
    S_padded = [S_windowed, zeros(1, N_fft - N_f)];

    % ---- IFFT to range domain ----
    range_profile = ifft(S_padded);
    range_profile_shifted = fftshift(range_profile);
    profile_dB = 20 * log10(max(abs(range_profile_shifted), 1e-30));

    % ---- Range axis ----
    B = max(freq) - min(freq);
    delta_r = c / (2 * B);
    r_max = c * N_f / (4 * B);
    range_axis = linspace(-r_max/2, r_max/2, N_fft);

    % ---- Plot ----
    plot(range_axis, profile_dB, 'b-', 'LineWidth', 1.8);
    xlabel('Range (m)', 'FontSize', 12);
    ylabel('|Range Profile| (dB)', 'FontSize', 12);
    title(sprintf(['Quick Range Profile (1D IFFT)  —  \\theta=%.1f°  \\phi=%.1f°  |  ' ...
                   'Window: %s  ZP: %dx'], ...
                   thVal, phVal, opt.windowFunc, opt.zeroPad), ...
          'FontSize', 10);

    % Mark peak scatterer
    [peakVal, peakI] = max(profile_dB);
    hold on;
    plot(range_axis(peakI), peakVal, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    text(range_axis(peakI), peakVal, ...
         sprintf('  %.1f dB @ %.3f m', peakVal, range_axis(peakI)), ...
         'FontSize', 9, 'Color', 'r');

    % Resolution indicator
    text(0.05, 0.92, sprintf('Range Res: %.3f m  |  BW: %.2f GHz  |  N_f: %d', ...
         delta_r, B/1e9, N_f), ...
         'Units', 'normalized', 'FontSize', 9, 'Color', [0.4 0.4 0.4]);

    grid on; box on;

    if opt.saveFigs
        saveFig(fig, fullfile(opt.saveDir, [opt.figPrefix '_quick_range_profile']), opt);
    end
end

%% ========================================================================
%  Figure 3: Spectrogram — |S| vs frequency vs angle
% ========================================================================
function plotSpectrogram(freq, fGHz, S, theta, phi, data, opt)
    fig = figure('Name', sprintf('FFE - Spectrogram [%s]', data.fileName), ...
                 'NumberTitle', 'off', ...
                 'Position', [50, 50, 1000, 450], 'Color', 'w');

    S_mag_dB = 20 * log10(max(abs(S), 1e-30));

    % Determine angle axis
    if data.it > data.ip
        % Sweep is primarily in theta
        angleAxis = theta(1, :);
        angLabel = '\theta (deg)';
    else
        angleAxis = phi(:, 1)';
        angLabel = '\phi (deg)';
    end

    imagesc(angleAxis, fGHz, S_mag_dB');
    set(gca, 'YDir', 'normal');
    xlabel(angLabel, 'FontSize', 12);
    ylabel('Frequency (GHz)', 'FontSize', 12);
    title(sprintf('Scattering Spectrogram  |S(\\theta, f)| (dB)  —  %s', data.fileName), ...
          'FontSize', 11);
    colormap(gca, 'jet');
    cb = colorbar;
    cb.Label.String = '|S| (dB)';
    cb.Label.FontSize = 11;
    axis tight;

    if opt.saveFigs
        saveFig(fig, fullfile(opt.saveDir, [opt.figPrefix '_spectrogram']), opt);
    end
end

%% ========================================================================
%  Figure 4: Frequency-Averaged RCS vs Angle
% ========================================================================
function plotFreqAveragedRCS(S, theta, phi, data, opt)
    fig = figure('Name', sprintf('FFE - Freq-Avg RCS [%s]', data.fileName), ...
                 'NumberTitle', 'off', ...
                 'Position', [80, 80, 800, 500], 'Color', 'w');

    % Average |S|^2 over frequency, then dB
    S_power = mean(abs(S).^2, 2);  % N_angles x 1
    S_power = reshape(S_power, size(theta));
    rcsAvg_dB = 10 * log10(max(S_power, 1e-30));

    if data.is2D
        % ---- 2D heatmap of freq-averaged RCS ----
        subplot(1, 2, 1);
        thetaDeg = theta(1, :);
        phiDeg   = phi(:, 1);
        [T, P] = meshgrid(thetaDeg, phiDeg);
        pcolor(T, P, rcsAvg_dB);
        shading interp;
        xlabel('\theta (deg)', 'FontSize', 12);
        ylabel('\phi (deg)', 'FontSize', 12);
        title('Frequency-Averaged RCS', 'FontSize', 11);
        colormap(gca, 'jet');
        colorbar;
        axis tight; grid on;
    end

    % ---- 1D center cut of averaged RCS ----
    subplot(1, 2, 2);
    if data.is2D
        midRow = round(data.ip / 2);
        angle = thetaDeg;
        rcsCut = rcsAvg_dB(midRow, :);
        angLabel = '\theta (deg)';
        sliceInfo = sprintf('\\phi = %.0f° (center cut)', phiDeg(midRow));
    else
        if data.ip == 1
            angle = theta(1, :);
            angLabel = '\theta (deg)';
        else
            angle = phi(:, 1)';
            angLabel = '\phi (deg)';
        end
        rcsCut = rcsAvg_dB(:)';
        sliceInfo = 'all angles';
    end

    plot(angle, rcsCut, 'b-', 'LineWidth', 2);
    xlabel(angLabel, 'FontSize', 12);
    ylabel('Freq-Averaged RCS (dBsm)', 'FontSize', 12);
    title(sprintf('Freq-Averaged RCS  [%s]', sliceInfo), 'FontSize', 11);
    grid on; box on;
    xlim([min(angle) max(angle)]);

    sgtitle(sprintf('Frequency-Averaged RCS Pattern  —  %s  (N_f = %d)', ...
            data.fileName, data.N_f), 'FontSize', 12);

    if opt.saveFigs
        saveFig(fig, fullfile(opt.saveDir, [opt.figPrefix '_freq_avg_rcs']), opt);
    end
end

%% ========================================================================
%  CREATEWINDOW  Generate a window function for IFFT pre-processing
% ========================================================================
function win = createWindow(N, windowType)
    switch lower(windowType)
        case 'hamming'
            win = hamming(N);
        case 'hann'
            win = hann(N);
        case 'blackman'
            win = blackman(N);
        case 'kaiser'
            win = kaiser(N, 2.5);
        case 'rect'
            win = ones(N, 1);
        otherwise
            warning('ffe_plot_frequency:UnknownWindow', ...
                    'Unknown window "%s". Using rectangular.', windowType);
            win = ones(N, 1);
    end
end

%% ========================================================================
function saveFig(fig, basePath, opt)
    if ~exist(opt.saveDir, 'dir')
        [status, msg] = mkdir(opt.saveDir);
        if ~status
            warning('ffe_plot_frequency:CannotCreateDir', ...
                    'Cannot create save directory: %s', msg);
            return;
        end
    end
    try
        exportgraphics(fig, [basePath '.png'], 'Resolution', 300);
        fprintf('  Figure saved: %s.png\n', basePath);
    catch
        try
            saveas(fig, [basePath '.png']);
            fprintf('  Figure saved: %s.png\n', basePath);
        catch ME
            warning('ffe_plot_frequency:SaveFailed', ...
                    'Failed to save figure: %s', ME.message);
        end
    end
end
