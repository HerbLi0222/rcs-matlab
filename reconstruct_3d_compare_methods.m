function reconstruct_3d_compare_methods(varargin)
% RECONSTRUCT_3D_COMPARE_METHODS  FFT vs PSF 两种 SAR 成像方法的三维重建对比
%
%   对同一个宽带散射数据集，在每个中心观测角度下分别使用：
%     - 方法一 (FFT)： 二维 IFFT 成像 + 峰值检测提取散射中心
%     - 方法二 (PSF)：  CLEAN 算法 + PSF 叠加提取散射中心
%
%   将两组散射中心分别反投影到三维空间，对比分析两种方法在
%   3D 空间中的重建差异。
%
%   方法原理:
%     - FFT 法: S(f,θ) → 2D IFFT → |I(r,x)| → 峰值检测 → (r,x,A)
%     - PSF 法: S(f,θ) → 2D IFFT → CLEAN 迭代 → PSF 叠加 → (r,x,A,φ)
%     - 3D 重建: P_3d = -r·k̂ - x·k̂_perp（修正相位约定符号翻转）
%
%   Usage:
%     >> reconstruct_3d_compare_methods
%     >> reconstruct_3d_compare_methods('DataFile', 'wideband_scattering_xxx.mat')
%     >> reconstruct_3d_compare_methods('CenterAngles', [30, 90, 150], ...
%           'SubApertureDeg', 10, 'FftThreshold', -6)

%% ========================================================================
%% 0. 参数解析
%% ========================================================================
addpath('lib');

p = inputParser;
p.addOptional('DataFile', '', @(x) ischar(x) || isstring(x) || isempty(x));
p.addParameter('CenterAngles', [30, 60, 90, 120, 150], @(x) isnumeric(x));
p.addParameter('SubApertureDeg', 10, @(x) isnumeric(x) && isscalar(x));
p.addParameter('WindowType', 'hamming', @(x) ischar(x) || isstring(x));
p.addParameter('ZeroPadFactor', 2, @(x) isnumeric(x) && isscalar(x));
p.addParameter('CleanMaxIter', 50, @(x) isnumeric(x) && isscalar(x));
p.addParameter('CleanThreshold', 0.05, @(x) isnumeric(x) && isscalar(x));
p.addParameter('CleanGain', 0.5, @(x) isnumeric(x) && isscalar(x));
p.addParameter('FftThreshold', -6, @(x) isnumeric(x) && isscalar(x));    % dB 阈值
p.addParameter('FftMinSep', 0.02, @(x) isnumeric(x) && isscalar(x));     % 最小峰间距 (m)
p.addParameter('FftMaxPeaks', 60, @(x) isnumeric(x) && isscalar(x));     % 每角度最大峰数
p.addParameter('MatchingRadius', 0.15, @(x) isnumeric(x) && isscalar(x)); % 匹配半径 (m)

if ~isempty(varargin) && (ischar(varargin{1}) || isstring(varargin{1})) ...
        && any(strcmp(char(varargin{1}), {'DataFile','CenterAngles',...
        'SubApertureDeg','WindowType','ZeroPadFactor','CleanMaxIter',...
        'CleanThreshold','CleanGain','FftThreshold','FftMinSep',...
        'FftMaxPeaks','MatchingRadius'}))
    varargin = [{''}, varargin];
end
p.parse(varargin{:});

resultFile        = p.Results.DataFile;
center_angles_deg = p.Results.CenterAngles;
subApertureDeg    = p.Results.SubApertureDeg;
windowType        = p.Results.WindowType;
zeroPadFactor     = p.Results.ZeroPadFactor;
cleanMaxIter      = p.Results.CleanMaxIter;
cleanThreshold    = p.Results.CleanThreshold;
cleanGain         = p.Results.CleanGain;
fftThreshold      = p.Results.FftThreshold;
fftMinSep         = p.Results.FftMinSep;
fftMaxPeaks       = p.Results.FftMaxPeaks;
matchingRadius    = p.Results.MatchingRadius;

%% ========================================================================
%% 1. 加载数据
%% ========================================================================
fprintf('========================================\n');
fprintf('  FFT vs PSF 3D Reconstruction Comparison\n');
fprintf('  (FFT 与 PSF 方法三维重建对比)\n');
fprintf('========================================\n\n');

% --- 1a. 加载 wideband_scattering 数据 ---
if isempty(resultFile)
    resultDir = fullfile('results');
    files = dir(fullfile(resultDir, 'wideband_scattering_*.mat'));
    if isempty(files)
        error('No wideband_scattering_*.mat found. Run main_wideband_scattering first.');
    end
    [~, idx] = sort([files.datenum], 'descend');
    resultFile = fullfile(resultDir, files(idx(1)).name);
end

fprintf('Loading: %s\n', resultFile);
data = load(resultFile);

S_complex   = data.S_complex;
freq_array  = data.freq_array(:)';
theta_array = data.theta_array;
phi_array   = data.phi_array;
c           = data.c;
B           = data.B;
delta_r     = data.delta_r;
ip          = data.ip;
it          = data.it;
N_f         = data.N_f;

f_center = (freq_array(1) + freq_array(end)) / 2;
lambda_c = c / f_center;

% 模型信息
if isfield(data, 'inputModel')
    [~, modelName, ~] = fileparts(data.inputModel);
    modelFile = fullfile('stl_models', data.inputModel);
else
    modelName = 'Unknown';
    modelFile = '';
end

if isfield(data, 'bbox_center')
    P_ref = data.bbox_center(:);
else
    P_ref = [0; 0; 0];
end

% 角度向量
if isfield(data, 'delt'), delt = data.delt; else delt = 1; end
if isfield(data, 'tstart'), tstart = data.tstart; else tstart = 0; end
if isfield(data, 'pstart'), pstart = data.pstart; else pstart = 0; end

if ip == 1 && it > 1
    theta_vec = theta_array(1, :)';
elseif it == 1 && ip > 1
    theta_vec = phi_array(:, 1);
else
    theta_vec = (1:size(S_complex, 1))';
end

fprintf('  Model: %s, fc=%.2f GHz, λc=%.4f m, B=%.1f GHz\n', ...
    modelName, f_center/1e9, lambda_c, B/1e9);
fprintf('  FFT threshold: %.0f dB,  CLEAN threshold: %.2f\n', ...
    fftThreshold, cleanThreshold);

% --- 1b. 加载 STL ---
[vertices, faces] = readSTLGeometry(modelFile);

%% ========================================================================
%% 2. 多角度双方法 SAR 成像
%% ========================================================================
N_angles = length(center_angles_deg);
fprintf('\n========================================\n');
fprintf('  Multi-Angle Dual-Method SAR Imaging\n');
fprintf('========================================\n');

% 存储两组散射中心
all_fft = [];   % [x, y, z, amplitude, theta_src]
all_psf = [];   % [x, y, z, amplitude, theta_src]

% 每角度统计
per_angle_stats = struct();

for i_ang = 1:N_angles
    theta0_deg = center_angles_deg(i_ang);
    half_span = subApertureDeg / 2;
    theta_min = theta0_deg - half_span;
    theta_max = theta0_deg + half_span;

    fprintf('\n--- Angle %d/%d: θ₀ = %.1f° [%.1f°, %.1f°] ---\n', ...
        i_ang, N_angles, theta0_deg, theta_min, theta_max);

    % --- 2a. 子孔径提取 ---
    angle_indices = find(theta_vec >= theta_min & theta_vec <= theta_max);
    if isempty(angle_indices)
        [~, sort_idx] = sort(abs(theta_vec - theta0_deg));
        n_take = max(2, floor(length(sort_idx) * subApertureDeg / ...
                               (theta_vec(end) - theta_vec(1))));
        angle_indices = sort(sort_idx(1:min(n_take, length(sort_idx))));
    end
    S_sub = S_complex(angle_indices, :);
    theta_sub = theta_vec(angle_indices);
    N_theta_sub = length(angle_indices);

    if N_theta_sub < 2
        warning('  Too few angles (%d), skipping.', N_theta_sub);
        continue;
    end

    % 实际角度跨度
    delta_theta_rad_actual = deg2rad(max(theta_sub) - min(theta_sub));
    delta_x = lambda_c / (2 * delta_theta_rad_actual);

    fprintf('  N_θ=%d, Δθ=%.1f° actual, Δx=%.3f m\n', ...
        N_theta_sub, rad2deg(delta_theta_rad_actual), delta_x);

    % --- 2b. 零填充 ---
    N_r_pad = N_f * zeroPadFactor;
    N_x_pad = N_theta_sub * zeroPadFactor;
    if mod(N_r_pad, 2) ~= 0, N_r_pad = N_r_pad + 1; end
    if mod(N_x_pad, 2) ~= 0, N_x_pad = N_x_pad + 1; end

    % --- 2c. 二维加窗 ---
    [W_2d, ~, ~] = create2DWindow(N_f, N_theta_sub, windowType);
    S_win = S_sub .* W_2d;

    % --- 2d. 补零 + 2D IFFT ---
    S_padded = padarray(S_win, [N_x_pad - N_theta_sub, N_r_pad - N_f], 0, 'post');
    I_fft = fftshift(ifft2(S_padded));

    % --- 2e. 坐标轴 ---
    delta_r_pad = delta_r * (N_f / N_r_pad);
    range_axis = (-floor(N_r_pad/2) : ceil(N_r_pad/2)-1) * delta_r_pad;
    delta_x_pad = delta_x * (N_theta_sub / N_x_pad);
    crossrange_axis = (-floor(N_x_pad/2) : ceil(N_x_pad/2)-1) * delta_x_pad;

    % --- 2f. PSF 函数 ---
    psf_func = @(r, x) sincf((2*B/c) * r) .* sincf((2*delta_theta_rad_actual/lambda_c) * x);

    %% ====================================================================
    %% 方法一：FFT 峰值检测
    %% ====================================================================
    I_fft_mag = abs(I_fft);
    I_fft_dB = 20 * log10(I_fft_mag / max(I_fft_mag(:)));

    [fft_peaks, ~] = detectPeaksFFT(I_fft_mag, I_fft_dB, range_axis, ...
        crossrange_axis, fftThreshold, fftMinSep, fftMaxPeaks);

    N_fft_peaks = size(fft_peaks, 1);
    fprintf('  Method 1 (FFT):  %d peaks detected\n', N_fft_peaks);

    %% ====================================================================
    %% 方法二：CLEAN + PSF
    %% ====================================================================
    [clean_components, ~] = cleanAlgorithm(I_fft, range_axis, crossrange_axis, ...
        psf_func, cleanThreshold, cleanMaxIter, cleanGain);

    % PSF 叠加合成图像
    [I_psf, ~] = buildPSFImage(clean_components, range_axis, crossrange_axis, psf_func);
    I_psf_dB = 20 * log10(abs(I_psf) / max(abs(I_psf(:))));

    N_psf_peaks = size(clean_components, 1);
    fprintf('  Method 2 (PSF): %d CLEAN centers\n', N_psf_peaks);

    %% ====================================================================
    %% 2D → 3D 坐标变换（修正相位约定）
    %% ====================================================================
    theta0_rad = deg2rad(theta0_deg);
    k_hat  = [sin(theta0_rad); 0; cos(theta0_rad)];
    k_perp = [cos(theta0_rad); 0; -sin(theta0_rad)];

    % --- FFT 峰 → 3D ---
    if ~isempty(fft_peaks)
        r_f = fft_peaks(:, 1);
        x_f = fft_peaks(:, 2);
        A_f = fft_peaks(:, 3);
        P_fft = - r_f * k_hat' - x_f * k_perp';
        all_fft = [all_fft; P_fft, A_f, theta0_deg * ones(N_fft_peaks, 1)]; %#ok<AGROW>
    end

    % --- CLEAN 中心 → 3D ---
    if ~isempty(clean_components)
        r_c = clean_components(:, 1);
        x_c = clean_components(:, 2);
        A_c = clean_components(:, 3);
        P_psf = - r_c * k_hat' - x_c * k_perp';
        all_psf = [all_psf; P_psf, A_c, theta0_deg * ones(N_psf_peaks, 1)]; %#ok<AGROW>
    end

    % 保存单角度统计
    per_angle_stats(i_ang).theta0 = theta0_deg;
    per_angle_stats(i_ang).N_fft = N_fft_peaks;
    per_angle_stats(i_ang).N_psf = N_psf_peaks;
    per_angle_stats(i_ang).I_fft_dB = I_fft_dB;
    per_angle_stats(i_ang).I_psf_dB = I_psf_dB;
    per_angle_stats(i_ang).range_axis = range_axis;
    per_angle_stats(i_ang).crossrange_axis = crossrange_axis;
    per_angle_stats(i_ang).fft_peaks = fft_peaks;
    per_angle_stats(i_ang).clean_components = clean_components;

end  % 角度循环

fprintf('\n========================================\n');
fprintf('  Summary: FFT total=%d, PSF total=%d\n', ...
    size(all_fft, 1), size(all_psf, 1));
fprintf('========================================\n');

%% ========================================================================
%% 3. 三维散射中心匹配与差异分析
%% ========================================================================
fprintf('\n--- 3D Matching Analysis ---\n');

[matched_pairs, dist_stats] = matchScatterers3D(all_fft, all_psf, matchingRadius);

fprintf('  Matching radius: %.3f m\n', matchingRadius);
fprintf('  FFT centers:        %d\n', size(all_fft, 1));
fprintf('  PSF centers:        %d\n', size(all_psf, 1));
fprintf('  Matched pairs:      %d\n', size(matched_pairs, 1));
fprintf('  FFT unmatched:      %d\n', dist_stats.N_fft_unmatched);
fprintf('  PSF unmatched:      %d\n', dist_stats.N_psf_unmatched);
if ~isempty(dist_stats.matched_dists)
    fprintf('  Match distances (m): min=%.4f, median=%.4f, mean=%.4f, max=%.4f\n', ...
        dist_stats.dist_min, dist_stats.dist_median, dist_stats.dist_mean, dist_stats.dist_max);
end

%% ========================================================================
%% 4. 可视化
%% ========================================================================
fprintf('\nGenerating comparison figures...\n');

% --- 通用绘图参数 ---
if ~isempty(all_fft)
    A_fft_norm = all_fft(:, 4) / max(all_fft(:, 4));
    sz_fft = 25 + 175 * A_fft_norm;
else
    sz_fft = [];
end
if ~isempty(all_psf)
    A_psf_norm = all_psf(:, 4) / max(all_psf(:, 4));
    sz_psf = 25 + 175 * A_psf_norm;
else
    sz_psf = [];
end

% 按来源角度着色
theta_fft = all_fft(:, 5);
theta_psf = all_psf(:, 5);
all_thetas = unique([center_angles_deg(:); theta_fft(:); theta_psf(:)]);
N_colors = length(center_angles_deg);
cmap = lines(N_colors);

color_fft = zeros(size(all_fft, 1), 3);
color_psf = zeros(size(all_psf, 1), 3);
for i = 1:size(all_fft, 1)
    [~, idx] = min(abs(center_angles_deg(:) - all_fft(i, 5)));
    color_fft(i, :) = cmap(idx, :);
end
for i = 1:size(all_psf, 1)
    [~, idx] = min(abs(center_angles_deg(:) - all_psf(i, 5)));
    color_psf(i, :) = cmap(idx, :);
end

%% --- 图1: 3D 并排对比 (FFT | PSF | Overlay) ---
figure(1);
clf;
set(gcf, 'Name', 'FFT vs PSF 3D Reconstruction Comparison', ...
         'NumberTitle', 'off', 'Position', [30, 60, 1400, 500]);

% --- 子图 1a: FFT 方法 ---
subplot(1, 3, 1);
hold on;
if ~isempty(vertices) && ~isempty(faces)
    patch('Faces', faces, 'Vertices', vertices, ...
          'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none', ...
          'FaceAlpha', 0.2);
end
if ~isempty(all_fft)
    for i = 1:size(all_fft, 1)
        scatter3(all_fft(i,1), all_fft(i,2), all_fft(i,3), ...
                 sz_fft(i), color_fft(i,:), 'filled', ...
                 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    end
end
scatter3(P_ref(1), P_ref(2), P_ref(3), 120, 'k', 'p', 'filled');
hold off;
axis equal; grid on; box on;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title(sprintf('Method 1: FFT + Peak Detection\n(%d centers)', size(all_fft,1)), ...
    'FontSize', 12, 'FontWeight', 'bold');
view(45, 25);
camlight('headlight'); lighting gouraud; material dull;

% --- 子图 1b: PSF/CLEAN 方法 ---
subplot(1, 3, 2);
hold on;
if ~isempty(vertices) && ~isempty(faces)
    patch('Faces', faces, 'Vertices', vertices, ...
          'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none', ...
          'FaceAlpha', 0.2);
end
if ~isempty(all_psf)
    for i = 1:size(all_psf, 1)
        scatter3(all_psf(i,1), all_psf(i,2), all_psf(i,3), ...
                 sz_psf(i), color_psf(i,:), 'filled', ...
                 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    end
end
scatter3(P_ref(1), P_ref(2), P_ref(3), 120, 'k', 'p', 'filled');
hold off;
axis equal; grid on; box on;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title(sprintf('Method 2: CLEAN + PSF\n(%d centers)', size(all_psf,1)), ...
    'FontSize', 12, 'FontWeight', 'bold');
view(45, 25);
camlight('headlight'); lighting gouraud; material dull;

% --- 子图 1c: 叠加对比 ---
subplot(1, 3, 3);
hold on;
if ~isempty(vertices) && ~isempty(faces)
    patch('Faces', faces, 'Vertices', vertices, ...
          'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none', ...
          'FaceAlpha', 0.2);
end
% FFT: 蓝色圆圈
if ~isempty(all_fft)
    scatter3(all_fft(:,1), all_fft(:,2), all_fft(:,3), ...
             max(20, sz_fft*0.6), [0.2 0.4 1.0], 'o', ...
             'LineWidth', 1.0, 'MarkerEdgeAlpha', 0.7);
end
% PSF: 红色十字
if ~isempty(all_psf)
    scatter3(all_psf(:,1), all_psf(:,2), all_psf(:,3), ...
             max(20, sz_psf*0.6), [1.0 0.2 0.2], '+', ...
             'LineWidth', 1.5, 'MarkerEdgeAlpha', 0.8);
end
% 匹配连线
if ~isempty(matched_pairs)
    for m = 1:size(matched_pairs, 1)
        p_fft = matched_pairs(m, 1:3);
        p_psf = matched_pairs(m, 4:6);
        plot3([p_fft(1) p_psf(1)], [p_fft(2) p_psf(2)], [p_fft(3) p_psf(3)], ...
              '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    end
end
scatter3(P_ref(1), P_ref(2), P_ref(3), 120, 'k', 'p', 'filled');
hold off;
axis equal; grid on; box on;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title(sprintf('Overlay: FFT (o) vs PSF (+)\n%d matched pairs', size(matched_pairs,1)), ...
    'FontSize', 12, 'FontWeight', 'bold');
legend({'FFT', 'PSF'}, 'Location', 'northeast', 'FontSize', 9);
view(45, 25);
camlight('headlight'); lighting gouraud; material dull;

sgtitle(sprintf('FFT vs PSF 3D Reconstruction | Model: %s', modelName), ...
       'FontSize', 14, 'FontWeight', 'bold');

%% --- 图2: 单角度 SAR 图像对比（选中间角度作为示例）---
figure(2);
clf;
set(gcf, 'Name', 'Per-Angle SAR Image Comparison', ...
         'NumberTitle', 'off', 'Position', [50, 80, 1300, 450]);

% 取中间角度展示
[~, mid_idx] = min(abs(center_angles_deg - median(center_angles_deg)));
N_show = min(3, N_angles);
show_indices = [];
if N_angles <= 3
    show_indices = 1:N_angles;
else
    show_indices = [1, mid_idx, N_angles];
end

for sp = 1:length(show_indices)
    s = show_indices(sp);
    if s > length(per_angle_stats), continue; end
    stat = per_angle_stats(s);

    % FFT 图像 + 检测到的峰
    subplot(2, length(show_indices), sp);
    imagesc(stat.range_axis, stat.crossrange_axis, stat.I_fft_dB);
    colormap('jet'); axis xy;
    hold on;
    if ~isempty(stat.fft_peaks)
        plot(stat.fft_peaks(:,1), stat.fft_peaks(:,2), 'wo', ...
             'MarkerSize', 6, 'LineWidth', 1.0);
    end
    hold off;
    xlabel('Range (m)'); ylabel('Cross-Range (m)');
    title(sprintf('FFT θ₀=%.0f° (%d peaks)', stat.theta0, stat.N_fft));
    caxis([-30, 0]); colorbar;

    % PSF 图像 + CLEAN 中心
    subplot(2, length(show_indices), length(show_indices) + sp);
    imagesc(stat.range_axis, stat.crossrange_axis, stat.I_psf_dB);
    colormap('jet'); axis xy;
    hold on;
    if ~isempty(stat.clean_components)
        amps = stat.clean_components(:,3);
        sz_c = 5 + 12 * amps / max(amps);
        for k = 1:size(stat.clean_components,1)
            plot(stat.clean_components(k,1), stat.clean_components(k,2), 'wo', ...
                 'MarkerSize', sz_c(k), 'LineWidth', 1.2);
        end
    end
    hold off;
    xlabel('Range (m)'); ylabel('Cross-Range (m)');
    title(sprintf('PSF+CLEAN θ₀=%.0f° (%d centers)', stat.theta0, stat.N_psf));
    caxis([-30, 0]); colorbar;
end

sgtitle(sprintf('SAR Image Comparison per Angle | Model: %s', modelName), ...
       'FontSize', 13, 'FontWeight', 'bold');

%% --- 图3: 统计分析 ---
figure(3);
clf;
set(gcf, 'Name', 'Statistical Analysis (统计分析)', ...
         'NumberTitle', 'off', 'Position', [100, 80, 1100, 700]);

% --- 3a: 每角度散射中心数量对比 ---
subplot(2, 3, 1);
theta_list = [per_angle_stats.theta0];
N_fft_list = [per_angle_stats.N_fft];
N_psf_list = [per_angle_stats.N_psf];
bar_width = 0.35;
b1 = bar(theta_list - bar_width/2, N_fft_list, bar_width, 'FaceColor', [0.2 0.4 1.0]);
hold on;
b2 = bar(theta_list + bar_width/2, N_psf_list, bar_width, 'FaceColor', [1.0 0.2 0.2]);
hold off;
xlabel('Center Angle θ₀ (deg)');
ylabel('Number of Centers');
title('Scatterer Count per Angle');
legend([b1, b2], {'FFT', 'PSF'}, 'Location', 'best');
grid on;

% --- 3b: 匹配距离直方图 ---
subplot(2, 3, 2);
if ~isempty(matched_pairs) && ~isempty(dist_stats.matched_dists)
    histogram(dist_stats.matched_dists, 25, 'FaceColor', [0.3 0.5 0.8], ...
              'EdgeColor', 'k', 'LineWidth', 0.5);
    hold on;
    xline(dist_stats.dist_median, 'r--', 'LineWidth', 2);
    xline(dist_stats.dist_mean, 'b--', 'LineWidth', 2);
    hold off;
    legend({sprintf('N=%d', length(dist_stats.matched_dists)), ...
           sprintf('Median=%.4f', dist_stats.dist_median), ...
           sprintf('Mean=%.4f', dist_stats.dist_mean)}, ...
           'Location', 'best', 'FontSize', 8);
end
xlabel('3D Distance (m)');
ylabel('Count');
title(sprintf('Match Distance Distribution (R=%.2f m)', matchingRadius));
grid on;

% --- 3c: X 坐标对比散点图 ---
subplot(2, 3, 3);
if ~isempty(matched_pairs)
    x_fft_m = matched_pairs(:, 1);
    x_psf_m = matched_pairs(:, 4);
    scatter(x_fft_m, x_psf_m, 20, [0.3 0.5 0.8], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    hold on;
    lims = [min([x_fft_m; x_psf_m]), max([x_fft_m; x_psf_m])];
    plot(lims, lims, 'k--', 'LineWidth', 1.5);
    hold off;
    xlabel('FFT X (m)'); ylabel('PSF X (m)');
    title('X-Coordinate: FFT vs PSF');
    axis equal; grid on;
    rho_x = corr(x_fft_m, x_psf_m);
    text(0.05, 0.95, sprintf('ρ=%.4f', rho_x), 'Units', 'normalized', ...
         'FontSize', 10, 'FontWeight', 'bold');
end

% --- 3d: Z 坐标对比散点图 ---
subplot(2, 3, 4);
if ~isempty(matched_pairs)
    z_fft_m = matched_pairs(:, 3);
    z_psf_m = matched_pairs(:, 6);
    scatter(z_fft_m, z_psf_m, 20, [0.8 0.3 0.3], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    hold on;
    lims = [min([z_fft_m; z_psf_m]), max([z_fft_m; z_psf_m])];
    plot(lims, lims, 'k--', 'LineWidth', 1.5);
    hold off;
    xlabel('FFT Z (m)'); ylabel('PSF Z (m)');
    title('Z-Coordinate: FFT vs PSF');
    axis equal; grid on;
    rho_z = corr(z_fft_m, z_psf_m);
    text(0.05, 0.95, sprintf('ρ=%.4f', rho_z), 'Units', 'normalized', ...
         'FontSize', 10, 'FontWeight', 'bold');
end

% --- 3e: 幅度对比 ---
subplot(2, 3, 5);
if ~isempty(matched_pairs)
    A_fft_m = matched_pairs(:, 7);
    A_psf_m = matched_pairs(:, 8);
    scatter(A_fft_m, A_psf_m, 20, [0.3 0.7 0.3], 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    hold on;
    lims = [0, max([A_fft_m; A_psf_m]) * 1.1];
    plot(lims, lims, 'k--', 'LineWidth', 1.5);
    hold off;
    xlabel('FFT Amplitude'); ylabel('PSF Amplitude');
    title('Amplitude: FFT vs PSF');
    axis equal; grid on;
    rho_a = corr(A_fft_m, A_psf_m);
    text(0.05, 0.95, sprintf('ρ=%.4f', rho_a), 'Units', 'normalized', ...
         'FontSize', 10, 'FontWeight', 'bold');
end

% --- 3f: 汇总表 ---
subplot(2, 3, 6);
axis off;
text_str = {
    sprintf('=== Comparison Summary ===');
    sprintf(' ');
    sprintf('FFT total centers:    %d', size(all_fft, 1));
    sprintf('PSF total centers:    %d', size(all_psf, 1));
    sprintf('Matched pairs:        %d', size(matched_pairs, 1));
    sprintf('Match rate (FFT):     %.1f%%', ...
        100 * size(matched_pairs, 1) / max(1, size(all_fft, 1)));
    sprintf('Match rate (PSF):     %.1f%%', ...
        100 * size(matched_pairs, 1) / max(1, size(all_psf, 1)));
    sprintf(' ');
    sprintf('Mean 3D distance:  %.4f m', dist_stats.dist_mean);
    sprintf('Median 3D distance: %.4f m', dist_stats.dist_median);
    sprintf('Max 3D distance:    %.4f m', dist_stats.dist_max);
    sprintf('(N matched pairs: %d)', size(matched_pairs, 1));
};
for t = 1:length(text_str)
    text(0.05, 0.95 - (t-1)*0.08, text_str{t}, 'Units', 'normalized', ...
         'FontSize', 10, 'FontName', 'FixedWidth');
end

sgtitle(sprintf('Statistical Analysis | Model: %s | %d Angles', ...
    modelName, N_angles), 'FontSize', 13, 'FontWeight', 'bold');

%% --- 图4: 三视图对比 (FFT 左列, PSF 右列) ---
figure(4);
clf;
set(gcf, 'Name', 'Three-View Comparison (三视图对比)', ...
         'NumberTitle', 'off', 'Position', [80, 60, 1100, 700]);

view_rows = {'XZ (Front)', 'YZ (Side)', 'XY (Top)'};
view_x = {1, 2, 1};  % 列索引
view_y = {3, 3, 2};  % 列索引
view_xl = {'X (m)', 'Y (m)', 'X (m)'};
view_yl = {'Z (m)', 'Z (m)', 'Y (m)'};

for v = 1:3
    % FFT
    subplot(3, 2, 2*v - 1);
    hold on;
    if ~isempty(vertices)
        bx = vertices(:, view_x{v}); by_p = vertices(:, view_y{v});
        try
            k_idx = boundary(bx, by_p, 0.5);
            if ~isempty(k_idx)
                fill(bx(k_idx), by_p(k_idx), [0.85 0.85 0.85], ...
                     'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 1);
            end
        catch
        end
    end
    if ~isempty(all_fft)
        for i = 1:size(all_fft, 1)
            scatter(all_fft(i, view_x{v}), all_fft(i, view_y{v}), ...
                    sz_fft(i)*0.5, color_fft(i,:), 'filled', ...
                    'MarkerEdgeColor', 'k', 'LineWidth', 0.2);
        end
    end
    scatter(P_ref(view_x{v}), P_ref(view_y{v}), 60, 'k', 'p', 'filled');
    hold off;
    axis equal; grid on;
    xlabel(view_xl{v}); ylabel(view_yl{v});
    title(sprintf('FFT: %s (%d centers)', view_rows{v}, size(all_fft,1)), ...
        'FontSize', 11, 'FontWeight', 'bold');

    % PSF
    subplot(3, 2, 2*v);
    hold on;
    if ~isempty(vertices)
        bx = vertices(:, view_x{v}); by_p = vertices(:, view_y{v});
        try
            k_idx = boundary(bx, by_p, 0.5);
            if ~isempty(k_idx)
                fill(bx(k_idx), by_p(k_idx), [0.85 0.85 0.85], ...
                     'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 1);
            end
        catch
        end
    end
    if ~isempty(all_psf)
        for i = 1:size(all_psf, 1)
            scatter(all_psf(i, view_x{v}), all_psf(i, view_y{v}), ...
                    sz_psf(i)*0.5, color_psf(i,:), 'filled', ...
                    'MarkerEdgeColor', 'k', 'LineWidth', 0.2);
        end
    end
    scatter(P_ref(view_x{v}), P_ref(view_y{v}), 60, 'k', 'p', 'filled');
    hold off;
    axis equal; grid on;
    xlabel(view_xl{v}); ylabel(view_yl{v});
    title(sprintf('PSF: %s (%d centers)', view_rows{v}, size(all_psf,1)), ...
        'FontSize', 11, 'FontWeight', 'bold');
end

sgtitle(sprintf('Three-View Comparison: FFT (left) vs PSF (right) | Model: %s', ...
    modelName), 'FontSize', 13, 'FontWeight', 'bold');

%% ========================================================================
%% 5. 保存结果
%% ========================================================================
fprintf('\nSaving results...\n');

nowStr = datestr(now, 'yyyymmddHHMMSS');

for figNum = 1:4
    if ishandle(figNum)
        figFile = fullfile('results', sprintf('compare3d_fig%d_%s.png', figNum, nowStr));
        saveas(figNum, figFile);
        fprintf('  Figure %d: %s\n', figNum, figFile);
    end
end

matFile = fullfile('results', ['compare3d_data_' nowStr '.mat']);
save(matFile, 'all_fft', 'all_psf', 'matched_pairs', 'dist_stats', ...
    'per_angle_stats', 'center_angles_deg', 'subApertureDeg', ...
    'fftThreshold', 'cleanThreshold', 'matchingRadius', ...
    'lambda_c', 'B', 'c', 'delta_r', 'modelName', '-v7.3');
fprintf('  Data: %s\n', matFile);

% 文本报告
rptFile = fullfile('results', ['compare3d_report_' nowStr '.txt']);
fid = fopen(rptFile, 'w');
fprintf(fid, '========================================\n');
fprintf(fid, 'FFT vs PSF 3D Reconstruction Comparison\n');
fprintf(fid, '========================================\n');
fprintf(fid, 'Model: %s\n', modelName);
fprintf(fid, 'Date:  %s\n', nowStr);
fprintf(fid, '\n--- Parameters ---\n');
fprintf(fid, 'Center angles: %s deg\n', num2str(center_angles_deg));
fprintf(fid, 'Sub-aperture: %.1f deg\n', subApertureDeg);
fprintf(fid, 'FFT threshold: %.0f dB\n', fftThreshold);
fprintf(fid, 'CLEAN threshold: %.3f, maxIter: %d, gain: %.2f\n', ...
    cleanThreshold, cleanMaxIter, cleanGain);
fprintf(fid, 'Matching radius: %.3f m\n', matchingRadius);
fprintf(fid, 'Phase correction: P_3d = -r*k_hat - x*k_perp\n');
fprintf(fid, '\n--- Per-Angle Counts ---\n');
fprintf(fid, '%-10s %-10s %-10s\n', 'Theta(deg)', 'N_FFT', 'N_PSF');
for i = 1:N_angles
    fprintf(fid, '%-10.1f %-10d %-10d\n', ...
        per_angle_stats(i).theta0, per_angle_stats(i).N_fft, per_angle_stats(i).N_psf);
end
fprintf(fid, '\n--- Summary ---\n');
fprintf(fid, 'FFT total:         %d\n', size(all_fft, 1));
fprintf(fid, 'PSF total:         %d\n', size(all_psf, 1));
fprintf(fid, 'Matched pairs:     %d\n', size(matched_pairs, 1));
if ~isempty(dist_stats.matched_dists)
    fprintf(fid, 'Mean distance:     %.4f m\n', dist_stats.dist_mean);
    fprintf(fid, 'Median distance:   %.4f m\n', dist_stats.dist_median);
    fprintf(fid, 'Max distance:      %.4f m\n', dist_stats.dist_max);
end
fclose(fid);
fprintf('  Report: %s\n', rptFile);

fprintf('\n========================================\n');
fprintf('  Comparison Complete\n');
fprintf('========================================\n');
fprintf('  FFT centers:        %d\n', size(all_fft, 1));
fprintf('  PSF centers:        %d\n', size(all_psf, 1));
fprintf('  Matched pairs:      %d\n', size(matched_pairs, 1));
if ~isempty(dist_stats.matched_dists)
    fprintf('  Mean 3D distance:   %.4f m\n', dist_stats.dist_mean);
    fprintf('  Median 3D distance: %.4f m\n', dist_stats.dist_median);
end
fprintf('========================================\n');

end  % reconstruct_3d_compare_methods

%% ========================================================================
%% 辅助函数
%% ========================================================================

% -------------------------------------------------------------------------
function [vertices, faces] = readSTLGeometry(modelFile)
    vertices = []; faces = [];
    if isempty(modelFile) || ~exist(modelFile, 'file')
        fprintf('  STL model not found: %s\n', modelFile);
        return;
    end
    try
        fid = fopen(modelFile, 'rb');
        if fid == -1, fprintf('  Cannot open STL: %s\n', modelFile); return; end
        fread(fid, 80, 'uint8=>uint8');
        nTri = fread(fid, 1, 'uint32=>uint32');
        allVerts = zeros(nTri * 3, 3);
        for i = 1:nTri
            fread(fid, 3, 'float32=>double');
            allVerts(3*i-2, :) = fread(fid, 3, 'float32=>double')';
            allVerts(3*i-1, :) = fread(fid, 3, 'float32=>double')';
            allVerts(3*i,   :) = fread(fid, 3, 'float32=>double')';
            fread(fid, 1, 'uint16=>uint16');
        end
        fclose(fid);
        [vertices, ~, ic] = unique(allVerts, 'rows', 'stable');
        faces = reshape(ic, 3, nTri)';
    catch ME
        fprintf('  STL read error: %s\n', ME.message);
    end
end

% -------------------------------------------------------------------------
function [W_2d, w_f, w_theta] = create2DWindow(N_f, N_theta, windowType)
    switch lower(windowType)
        case 'hamming'
            w_f = hamming(N_f)'; w_theta = hamming(N_theta);
        case 'hanning'
            w_f = hanning(N_f)'; w_theta = hanning(N_theta);
        case 'blackman'
            w_f = blackman(N_f)'; w_theta = blackman(N_theta);
        case 'kaiser'
            w_f = kaiser(N_f, 2.5)'; w_theta = kaiser(N_theta, 2.5);
        case 'none'
            w_f = ones(1, N_f); w_theta = ones(N_theta, 1);
        otherwise
            w_f = hamming(N_f)'; w_theta = hamming(N_theta);
    end
    w_f = w_f(:)';
    w_theta = w_theta(:);
    W_2d = w_theta * w_f;
end

% -------------------------------------------------------------------------
function [peaks, I_smoothed] = detectPeaksFFT(I_mag, I_dB, range_axis, ...
    crossrange_axis, threshold_dB, min_sep, max_peaks)
% 从 FFT 幅度图像中检测散射峰
%
% Input:
%   I_mag          - FFT 幅度图像 (N_x × N_r)
%   I_dB           - dB 归一化图像
%   range_axis     - 距离轴
%   crossrange_axis - 方位轴
%   threshold_dB   - dB 阈值（相对于 0dB 峰值）
%   min_sep        - 最小峰间距 (m)
%   max_peaks      - 最大峰数
% Output:
%   peaks          - N × 3: [range, cross_range, amplitude]

    peaks = zeros(0, 3);

    N_r = length(range_axis);
    N_x = length(crossrange_axis);

    % 高斯平滑
    I_smoothed = imgaussfilt(I_dB, 0.6);

    % 找局部最大值（3×3邻域）
    local_max = zeros(size(I_smoothed));
    for i = 2:N_x-1
        for j = 2:N_r-1
            patch = I_smoothed(i-1:i+1, j-1:j+1);
            if I_smoothed(i,j) == max(patch(:)) && I_smoothed(i,j) > threshold_dB
                local_max(i,j) = 1;
            end
        end
    end

    [peak_rows, peak_cols] = find(local_max);
    N_raw = length(peak_rows);

    if N_raw == 0
        % 放宽阈值重试
        relaxed_thresh = threshold_dB - 6;
        for i = 2:N_x-1
            for j = 2:N_r-1
                patch = I_smoothed(i-1:i+1, j-1:j+1);
                if I_smoothed(i,j) == max(patch(:)) && I_smoothed(i,j) > relaxed_thresh
                    local_max(i,j) = 1;
                end
            end
        end
        [peak_rows, peak_cols] = find(local_max);
        N_raw = length(peak_rows);
    end

    if N_raw == 0, return; end

    % 提取峰值属性 (含亚像素精化)
    raw_peaks = zeros(N_raw, 3);
    for k = 1:N_raw
        r_idx = peak_cols(k);
        x_idx = peak_rows(k);

        % 二次插值亚像素精化
        if r_idx > 1 && r_idx < N_r
            v_l = I_mag(x_idx, r_idx-1);
            v_c = I_mag(x_idx, r_idx);
            v_r = I_mag(x_idx, r_idx+1);
            delta_r_sub = (v_l - v_r) / (2*(v_l - 2*v_c + v_r + 1e-15));
            r_fine = range_axis(r_idx) + delta_r_sub * (range_axis(2) - range_axis(1));
        else
            r_fine = range_axis(r_idx);
        end

        if x_idx > 1 && x_idx < N_x
            v_d = I_mag(x_idx-1, r_idx);
            v_c = I_mag(x_idx, r_idx);
            v_u = I_mag(x_idx+1, r_idx);
            delta_x_sub = (v_d - v_u) / (2*(v_d - 2*v_c + v_u + 1e-15));
            x_fine = crossrange_axis(x_idx) + delta_x_sub * (crossrange_axis(2) - crossrange_axis(1));
        else
            x_fine = crossrange_axis(x_idx);
        end

        raw_peaks(k, :) = [r_fine, x_fine, I_mag(x_idx, r_idx)];
    end

    % 按幅度降序排列
    raw_peaks = sortrows(raw_peaks, -3);

    % 合并过于靠近的峰
    merged = raw_peaks(1, :);
    for k = 2:N_raw
        too_close = false;
        for m = 1:size(merged, 1)
            dr = raw_peaks(k, 1) - merged(m, 1);
            dx = raw_peaks(k, 2) - merged(m, 2);
            if sqrt(dr^2 + dx^2) < min_sep
                too_close = true;
                break;
            end
        end
        if ~too_close
            merged = [merged; raw_peaks(k, :)]; %#ok<AGROW>
        end
        if size(merged, 1) >= max_peaks, break; end
    end

    peaks = merged;
end

% -------------------------------------------------------------------------
function [components, residual] = cleanAlgorithm(I_complex, range_axis, ...
    crossrange_axis, psf_func, threshold, maxIter, gain)
% CLEAN 迭代散射中心提取

    residual = I_complex;
    residual_mag = abs(residual);
    components = [];
    initial_peak = max(residual_mag(:));
    stop_threshold = threshold * initial_peak;

    N_r = length(range_axis);
    N_x = length(crossrange_axis);

    for iter = 1:maxIter
        [peak_val, peak_idx] = max(residual_mag(:));
        if peak_val < stop_threshold, break; end

        [peak_row, peak_col] = ind2sub([N_x, N_r], peak_idx);
        r_peak = range_axis(peak_col);
        x_peak = crossrange_axis(peak_row);
        complex_peak = residual(peak_row, peak_col);
        A_k = abs(complex_peak);
        phi_k = angle(complex_peak);

        components = [components; r_peak, x_peak, A_k, phi_k]; %#ok<AGROW>

        [R_mesh, X_mesh] = meshgrid(range_axis, crossrange_axis);
        psf_shifted = psf_func(R_mesh - r_peak, X_mesh - x_peak);
        residual = residual - gain * complex_peak * psf_shifted;
        residual_mag = abs(residual);
    end
end

% -------------------------------------------------------------------------
function [I_psf, I_psf_mag] = buildPSFImage(components, range_axis, ...
    crossrange_axis, psf_func)
% PSF 叠加合成图像

    N_r = length(range_axis);
    N_x = length(crossrange_axis);
    I_psf = zeros(N_x, N_r);

    if isempty(components)
        I_psf_mag = abs(I_psf);
        return;
    end

    for k = 1:size(components, 1)
        r_k = components(k, 1);
        x_k = components(k, 2);
        A_k = components(k, 3);
        phi_k = components(k, 4);
        gain_c = A_k * exp(1j * phi_k);

        [R_mesh, X_mesh] = meshgrid(range_axis, crossrange_axis);
        psf_shifted = psf_func(R_mesh - r_k, X_mesh - x_k);
        I_psf = I_psf + gain_c * psf_shifted;
    end
    I_psf_mag = abs(I_psf);
end

% -------------------------------------------------------------------------
function [matched_pairs, stats] = matchScatterers3D(all_fft, all_psf, radius)
% 在 3D 空间中匹配 FFT 和 PSF 散射中心（最近邻双向匹配）
%
% Output:
%   matched_pairs - M × 8: [fft_x, fft_y, fft_z, psf_x, psf_y, psf_z, fft_A, psf_A]
%   stats         - 包含匹配距离统计的结构体

    matched_pairs = zeros(0, 8);
    stats = struct('matched_dists', [], 'N_fft_unmatched', 0, 'N_psf_unmatched', 0, ...
                   'dist_min', NaN, 'dist_median', NaN, 'dist_mean', NaN, 'dist_max', NaN);

    N_fft = size(all_fft, 1);
    N_psf = size(all_psf, 1);

    if N_fft == 0 || N_psf == 0
        stats.N_fft_unmatched = N_fft;
        stats.N_psf_unmatched = N_psf;
        return;
    end

    % 距离矩阵
    pos_fft = all_fft(:, 1:3);
    pos_psf = all_psf(:, 1:3);
    D = pdist2(pos_fft, pos_psf);

    % 双向最近邻匹配
    used_fft = false(N_fft, 1);
    used_psf = false(N_psf, 1);
    dists = [];

    while true
        % 找全局最小距离
        [min_val, min_idx] = min(D(:));
        if min_val > radius, break; end

        [i_fft, i_psf] = ind2sub(size(D), min_idx);
        if used_fft(i_fft) || used_psf(i_psf)
            D(i_fft, i_psf) = inf;  % 已用，跳过
            continue;
        end

        % 确认双向最近邻
        [~, nn_psf] = min(D(i_fft, :));
        [~, nn_fft] = min(D(:, i_psf));
        if nn_psf == i_psf && nn_fft == i_fft
            % 双向匹配成功
            matched_pairs = [matched_pairs; ...
                all_fft(i_fft, 1:3), all_psf(i_psf, 1:3), ...
                all_fft(i_fft, 4), all_psf(i_psf, 4)]; %#ok<AGROW>
            dists = [dists; min_val]; %#ok<AGROW>
            used_fft(i_fft) = true;
            used_psf(i_psf) = true;
            D(i_fft, :) = inf;
            D(:, i_psf) = inf;
        else
            D(i_fft, i_psf) = inf;
        end
    end

    stats.matched_dists = dists;
    stats.N_fft_unmatched = sum(~used_fft);
    stats.N_psf_unmatched = sum(~used_psf);
    if ~isempty(dists)
        stats.dist_min = min(dists);
        stats.dist_median = median(dists);
        stats.dist_mean = mean(dists);
        stats.dist_max = max(dists);
    end
end

% -------------------------------------------------------------------------
function y = sincf(x)
    y = ones(size(x));
    nonzero = abs(x) > 1e-15;
    y(nonzero) = sin(pi * x(nonzero)) ./ (pi * x(nonzero));
end
