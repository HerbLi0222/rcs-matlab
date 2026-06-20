function main_sar_imaging(varargin)
% MAIN_SAR_IMAGING  合成孔径成像 (Synthetic Aperture Radar Imaging)
%
%   基于 wideband_scattering 数据 S(f, theta) 进行合成孔径成像。
%   实现两种方法：
%     - 方法一：传统二维FFT成像（小角度近似）
%     - 方法二：PSF（点扩散函数）叠加成像（CLEAN算法提取散射中心）
%
%   实验内容：
%     1. 数据加载与子孔径裁剪（全孔径、1/2孔径、1/4孔径）
%     2. 二维加窗预处理
%     3. FFT法成像（直接2D IFFT）
%     4. PSF叠加法成像（CLEAN → PSF构建 → 叠加合成）
%     5. 多孔径对比分析（CC、SSIM、展宽宽度）
%
%   输入文件:
%     results/wideband_scattering_<timestamp>.mat  (由 main_wideband_scattering 生成)
%
%   参考文档:
%     合成孔径成像方法.txt
%
%   Usage:
%     >> main_sar_imaging                % 自动查找最新 wideband_scattering 数据
%     >> main_sar_imaging('file.mat')    % 指定数据文件
%     >> main_sar_imaging(..., 'CenterAngle', 90)  % 指定子孔径中心角 (默认: 90°)

%% ========================================================================
%% 0. 参数解析
%% ========================================================================
addpath('lib');

% 参数名列表（用于区分"文件路径"和"名称-值对参数"）
paramNames = {'DataFile','CenterAngle','AperturesDeg','WindowType',...
    'ZeroPadFactor','CleanThreshold','CleanMaxIter','CleanGain'};

% 默认参数：DataFile 为可选位置参数（第一个），其余为名称-值对
p = inputParser;
p.addOptional('DataFile', '', @(x) ischar(x) || isstring(x) || isempty(x));
p.addParameter('CenterAngle', 90, @(x) isnumeric(x) && isscalar(x));
p.addParameter('AperturesDeg', [10, 5, 2.5], @(x) isnumeric(x));
p.addParameter('WindowType', 'hamming', @(x) ischar(x) || isstring(x));
p.addParameter('ZeroPadFactor', 2, @(x) isnumeric(x) && isscalar(x));
p.addParameter('CleanThreshold', 0.05, @(x) isnumeric(x) && isscalar(x));
p.addParameter('CleanMaxIter', 50, @(x) isnumeric(x) && isscalar(x));
p.addParameter('CleanGain', 0.5, @(x) isnumeric(x) && isscalar(x));

% 如果第一个位置参数恰好是某个参数名 → 当作名称-值对处理，DataFile 用默认值
if ~isempty(varargin) && (ischar(varargin{1}) || isstring(varargin{1})) ...
        && any(strcmp(char(varargin{1}), paramNames))
    varargin = [{''}, varargin];  % 空 DataFile 占位
end
p.parse(varargin{:});

resultFile       = p.Results.DataFile;
theta_center_deg = p.Results.CenterAngle;
apertures_deg    = p.Results.AperturesDeg;
windowType       = p.Results.WindowType;
zeroPadFactor    = p.Results.ZeroPadFactor;
cleanThreshold   = p.Results.CleanThreshold;
cleanMaxIter     = p.Results.CleanMaxIter;
cleanGain        = p.Results.CleanGain;

%% ========================================================================
%% 1. 数据加载
%% ========================================================================
fprintf('========================================\n');
fprintf('  Synthetic Aperture Radar Imaging\n');
fprintf('  (合成孔径成像)\n');
fprintf('========================================\n\n');

% 自动查找最新 wideband_scattering 文件
if isempty(resultFile)
    resultFile = findLatestResultFile('wideband_scattering_*.mat');
    if isempty(resultFile)
        error('main_sar_imaging:NoDataFile', ...
              'No wideband_scattering_*.mat found in results/.\nRun main_wideband_scattering first.');
    end
end

fprintf('Loading wideband scattering data...\n');
fprintf('  File: %s\n', resultFile);

if ~exist(resultFile, 'file')
    error('main_sar_imaging:FileNotFound', 'File not found: %s', resultFile);
end

data = load(resultFile);

% 验证必要字段
requiredFields = {'S_complex', 'freq_array', 'theta_array', 'phi_array', ...
                  'c', 'B', 'delta_r', 'ip', 'it', 'N_f', 'N_angles_total'};
missingFields = {};
for k = 1:length(requiredFields)
    if ~isfield(data, requiredFields{k})
        missingFields{end+1} = requiredFields{k}; %#ok<AGROW>
    end
end
if ~isempty(missingFields)
    error('main_sar_imaging:MissingField', ...
          'Missing required field(s): %s.', strjoin(missingFields, ', '));
end

% 提取数据
S_complex      = data.S_complex;       % 复散射场矩阵 (N_angles × N_f)
freq_array     = data.freq_array(:)';  % 频率数组 (Hz), 行向量
theta_array    = data.theta_array;     % theta 角度数组
phi_array      = data.phi_array;       % phi 角度数组
c              = data.c;               % 光速 (m/s)
B              = data.B;               % 带宽 (Hz)
delta_r        = data.delta_r;         % 距离分辨率 (m)
ip             = data.ip;              % phi 步数
it             = data.it;              % theta 步数
N_f            = data.N_f;             % 频点数
N_angles_total = data.N_angles_total;  % 总角度数

% 模型名
if isfield(data, 'inputModel')
    [~, modelName, ~] = fileparts(data.inputModel);
else
    modelName = 'Unknown';
end

% 角度步长
if isfield(data, 'delt'), delt = data.delt; else delt = 1; end
if isfield(data, 'delp'), delp = data.delp; else delp = 1; end
if isfield(data, 'tstart'), tstart = data.tstart; else tstart = 0; end
if isfield(data, 'pstart'), pstart = data.pstart; else pstart = 0; end

% 计算中心频率和波长
f_center = (freq_array(1) + freq_array(end)) / 2;
lambda_c = c / f_center;

fprintf('  Model:          %s\n', modelName);
fprintf('  Frequency:      %.2f - %.2f GHz\n', freq_array(1)/1e9, freq_array(end)/1e9);
fprintf('  Bandwidth B:    %.2f GHz\n', B/1e9);
fprintf('  Center freq fc: %.2f GHz\n', f_center/1e9);
fprintf('  Wavelength λc:  %.4f m\n', lambda_c);
fprintf('  Range res Δr:   %.4f m\n', delta_r);
fprintf('  Matrix size:    %d angles × %d frequencies\n', N_angles_total, N_f);
fprintf('  Theta scan:     %.1f:%.1f:%.1f deg\n', tstart, delt, tstart + (it-1)*delt);

% 构建 theta 角度向量
if ip == 1 && it > 1
    % 纯 theta 扫描
    theta_vec = theta_array(1, :)';
    angleLabel = '\theta';
elseif it == 1 && ip > 1
    % 纯 phi 扫描
    theta_vec = phi_array(:, 1);
    angleLabel = '\phi';
else
    % 二维扫描
    theta_vec = (1:N_angles_total)';
    angleLabel = 'idx';
end

fprintf('  Angle vector:   %d points, %.1f° to %.1f°\n', ...
    length(theta_vec), theta_vec(1), theta_vec(end));

%% ========================================================================
%% 2. 成像参数配置
%% ========================================================================
fprintf('\n========================================\n');
fprintf('  SAR Imaging Configuration\n');
fprintf('========================================\n');

% 子孔径配置
N_apertures = length(apertures_deg);
aperture_labels = cell(1, N_apertures);
for i = 1:N_apertures
    aperture_labels{i} = sprintf('%.1f°', apertures_deg(i));
end

fprintf('  Center angle θ0: %.1f°\n', theta_center_deg);
fprintf('  Sub-apertures:   %s\n', strjoin(aperture_labels, ', '));
fprintf('  Window:          %s\n', windowType);
fprintf('  Zero-pad factor: %dx\n', zeroPadFactor);
fprintf('  CLEAN threshold: %.2f\n', cleanThreshold);
fprintf('  CLEAN max iter:  %d\n', cleanMaxIter);
fprintf('  CLEAN gain:      %.2f\n', cleanGain);

%% ========================================================================
%% 3. 对每个子孔径进行 SAR 成像
%% ========================================================================
fprintf('\n========================================\n');
fprintf('  SAR Imaging per Sub-Aperture\n');
fprintf('========================================\n');

% 存储各孔径的结果
results_fft    = cell(1, N_apertures);
results_psf    = cell(1, N_apertures);
results_clean  = cell(1, N_apertures);
aperture_info  = cell(1, N_apertures);

for i_ap = 1:N_apertures
    delta_theta_deg = apertures_deg(i_ap);
    half_span = delta_theta_deg / 2;
    theta_min = theta_center_deg - half_span;
    theta_max = theta_center_deg + half_span;

    fprintf('\n--- Sub-Aperture %d/%d: Δθ = %.1f° [%.1f°, %.1f°] ---\n', ...
        i_ap, N_apertures, delta_theta_deg, theta_min, theta_max);

    %% 3a. 提取子孔径数据
    [S_sub, theta_sub, N_theta_sub] = extractSubAperture(...
        S_complex, theta_vec, theta_min, theta_max);

    if N_theta_sub < 2
        warning('  Too few angles (%d), skipping this aperture.', N_theta_sub);
        continue;
    end

    % 转换为弧度
    delta_theta_rad = deg2rad(delta_theta_deg);

    % 方位分辨率
    delta_x = lambda_c / (2 * delta_theta_rad);

    fprintf('  N_theta = %d angles, Δθ_rad = %.4f rad\n', N_theta_sub, delta_theta_rad);
    fprintf('  Cross-range resolution Δx = λc/(2Δθ) = %.4f m\n', delta_x);

    % 零填充
    N_r_pad = N_f * zeroPadFactor;
    N_x_pad = N_theta_sub * zeroPadFactor;
    if mod(N_r_pad, 2) ~= 0, N_r_pad = N_r_pad + 1; end
    if mod(N_x_pad, 2) ~= 0, N_x_pad = N_x_pad + 1; end

    %% 3b. 二维加窗
    [W_2d, ~, ~] = create2DWindow(N_f, N_theta_sub, windowType);

    S_win = S_sub .* W_2d;  % N_theta_sub × N_f

    %% 3c. 坐标轴标定
    % 距离向 (Range)
    delta_r_pad = delta_r * (N_f / N_r_pad);
    range_axis = (-floor(N_r_pad/2) : ceil(N_r_pad/2)-1) * delta_r_pad;

    % 方位向 (Cross-Range)
    delta_x_pad = delta_x * (N_theta_sub / N_x_pad);
    crossrange_axis = (-floor(N_x_pad/2) : ceil(N_x_pad/2)-1) * delta_x_pad;

    %% ====================================================================
    %% 方法一：二维FFT成像
    %% ====================================================================
    fprintf('  Method 1: 2D FFT imaging...\n');
    tic;

    % 补零到成像网格
    S_padded = padarray(S_win, [N_x_pad - N_theta_sub, N_r_pad - N_f], 0, 'post');

    % 二维 IFFT
    I_fft_complex = fftshift(ifft2(S_padded));
    I_fft_mag = abs(I_fft_complex);

    % 归一化 (0 dB)
    I_fft_dB = 20 * log10(I_fft_mag / max(I_fft_mag(:)));

    t_fft = toc;
    fprintf('    Done in %.3f s\n', t_fft);

    %% ====================================================================
    %% 方法二：PSF叠加成像（CLEAN算法）
    %% ====================================================================
    fprintf('  Method 2: PSF superposition (CLEAN)...\n');
    tic;

    % PSF解析式: h(r, x) = sinc(2B/c * r) * sinc(2*Δθ/λ_c * x)
    % 其中 sinc(u) = sin(pi*u) / (pi*u)
    psf_func = @(r, x) sincf((2*B/c) * r) .* sincf((2*delta_theta_rad/lambda_c) * x);

    % CLEAN 算法提取散射中心
    [clean_components, ~] = cleanAlgorithm(...
        I_fft_complex, range_axis, crossrange_axis, psf_func, ...
        cleanThreshold, cleanMaxIter, cleanGain);

    N_clean = size(clean_components, 1);
    fprintf('    Extracted %d scattering centers\n', N_clean);

    % PSF 叠加合成图像
    [I_psf_complex, I_psf_mag] = buildPSFImage(...
        clean_components, range_axis, crossrange_axis, psf_func);

    % 归一化 (0 dB)
    I_psf_dB = 20 * log10(I_psf_mag / max(I_psf_mag(:)));

    t_psf = toc;
    fprintf('    Done in %.3f s\n', t_psf);

    %% 3d. 保存结果
    results_fft{i_ap} = struct(...
        'dB',         I_fft_dB, ...
        'complex',    I_fft_complex, ...
        'mag',        I_fft_mag, ...
        'range_axis', range_axis, ...
        'crossrange_axis', crossrange_axis);

    results_psf{i_ap} = struct(...
        'dB',         I_psf_dB, ...
        'complex',    I_psf_complex, ...
        'mag',        I_psf_mag, ...
        'range_axis', range_axis, ...
        'crossrange_axis', crossrange_axis);

    results_clean{i_ap} = clean_components;

    aperture_info{i_ap} = struct(...
        'delta_theta_deg', delta_theta_deg, ...
        'delta_theta_rad', delta_theta_rad, ...
        'delta_x',         delta_x, ...
        'N_theta_sub',     N_theta_sub, ...
        'N_r_pad',         N_r_pad, ...
        'N_x_pad',         N_x_pad, ...
        'delta_r_pad',     delta_r_pad, ...
        'delta_x_pad',     delta_x_pad, ...
        'theta_min',       theta_min, ...
        'theta_max',       theta_max);

end  % 孔径循环

%% ========================================================================
%% 4. 多孔径对比分析
%% ========================================================================
fprintf('\n========================================\n');
fprintf('  Multi-Aperture Comparison Analysis\n');
fprintf('========================================\n');

% 去除跳过的孔径
valid_idx = find(~cellfun(@isempty, results_fft));
N_valid = length(valid_idx);

if N_valid < 2
    warning('  Need at least 2 valid apertures for comparison.');
end

% 计算对比指标 + 逐方法指标
metrics = struct();
for i = 1:N_valid
    idx = valid_idx(i);

    % --- 对比指标 ---
    cc_val   = computeCorrelationCoefficient(results_fft{idx}.mag, results_psf{idx}.mag);
    ssim_val = computeSSIM(results_fft{idx}.dB, results_psf{idx}.dB);

    % --- 逐方法: 展宽宽度 ---
    bw_fft = measureBeamwidth(results_fft{idx}.dB);
    bw_psf = measureBeamwidth(results_psf{idx}.dB);

    % --- 逐方法: 散射点数（FFT: -20dB 以上像素; PSF: CLEAN 中心数）---
    n_fft = sum(results_fft{idx}.dB(:) > -20);
    n_psf = size(results_clean{idx}, 1);

    % --- 逐方法: 平均幅度 ---
    mask_fft = results_fft{idx}.dB > -20;
    if any(mask_fft)
        mean_amp_fft = mean(results_fft{idx}.mag(mask_fft));
    else
        mean_amp_fft = 0;
    end
    if ~isempty(results_clean{idx})
        mean_amp_psf = mean(results_clean{idx}(:, 3));
    else
        mean_amp_psf = 0;
    end

    metrics(i).cc       = cc_val;
    metrics(i).ssim     = ssim_val;
    metrics(i).bw_fft   = bw_fft;
    metrics(i).bw_psf   = bw_psf;
    metrics(i).n_fft    = n_fft;
    metrics(i).n_psf    = n_psf;
    metrics(i).amp_fft  = mean_amp_fft;
    metrics(i).amp_psf  = mean_amp_psf;
    metrics(i).delta_theta = aperture_info{idx}.delta_theta_deg;
    metrics(i).delta_x     = aperture_info{idx}.delta_x;

    fprintf('  Δθ=%.1f°: CC=%.3f SSIM=%.3f | BW FFT=%.1f/PSF=%.1f | N FFT=%d/PSF=%d | Amp FFT=%.2e/PSF=%.2e\n', ...
        metrics(i).delta_theta, cc_val, ssim_val, bw_fft, bw_psf, n_fft, n_psf, mean_amp_fft, mean_amp_psf);
end

%% ========================================================================
%% 5. 可视化
%% ========================================================================
fprintf('\nGenerating figures...\n');

% 全局色彩范围 (跨所有孔径)
all_dB_vals = [];
for i = 1:N_valid
    idx = valid_idx(i);
    all_dB_vals = [all_dB_vals; results_fft{idx}.dB(:); results_psf{idx}.dB(:)]; %#ok<AGROW>
end
caxis_min = prctile(all_dB_vals, 1);
caxis_max = max(all_dB_vals(:));
caxis_range = [caxis_min, caxis_max];

%% --- 图1: FFT 成像结果（多孔径并列）---
figure(1);
clf;
set(gcf, 'Name', 'FFT SAR Imaging Results (二维FFT合成孔径成像)', ...
         'NumberTitle', 'off', 'Position', [50, 80, 300*N_valid, 500]);

for i = 1:N_valid
    idx = valid_idx(i);
    subplot(2, N_valid, i);

    r_ax = results_fft{idx}.range_axis;
    x_ax = results_fft{idx}.crossrange_axis;
    imagesc(r_ax, x_ax, results_fft{idx}.dB);
    colormap('jet');
    axis xy;
    xlabel('Range (m)');
    ylabel('Cross-Range (m)');
    title(sprintf('FFT Δθ=%.1f°', aperture_info{idx}.delta_theta_deg));
    caxis(caxis_range);
    colorbar;

    subplot(2, N_valid, N_valid + i);
    imagesc(r_ax, x_ax, results_fft{idx}.dB);
    colormap('jet');
    axis xy;
    xlabel('Range (m)');
    ylabel('Cross-Range (m)');
    title(sprintf('FFT Δθ=%.1f° (zoomed)', aperture_info{idx}.delta_theta_deg));
    caxis(caxis_range);
    colorbar;
    % 缩放到目标区域
    xlim([r_ax(round(end*0.3)), r_ax(round(end*0.7))]);
end

sgtitle(sprintf('Method 1: 2D FFT SAR Imaging | Model: %s | θ_0=%.0f°', ...
    modelName, theta_center_deg), 'FontSize', 13, 'FontWeight', 'bold');

%% --- 图2: PSF 成像结果（多孔径并列）---
figure(2);
clf;
set(gcf, 'Name', 'PSF SAR Imaging Results (PSF叠加合成孔径成像)', ...
         'NumberTitle', 'off', 'Position', [100, 80, 300*N_valid, 500]);

for i = 1:N_valid
    idx = valid_idx(i);
    r_ax = results_psf{idx}.range_axis;
    x_ax = results_psf{idx}.crossrange_axis;
    subplot(2, N_valid, i);
    imagesc(r_ax, x_ax, results_psf{idx}.dB);
    colormap('jet');
    axis xy;
    xlabel('Range (m)');
    ylabel('Cross-Range (m)');
    title(sprintf('PSF Δθ=%.1f°', aperture_info{idx}.delta_theta_deg));
    caxis(caxis_range);
    colorbar;

    % 标注 CLEAN 散射中心位置
    subplot(2, N_valid, N_valid + i);
    hold on;
    imagesc(r_ax, x_ax, results_psf{idx}.dB);
    colormap('jet');
    axis xy;
    comps = results_clean{idx};
    if ~isempty(comps)
        for k = 1:size(comps, 1)
            plot(comps(k, 1), comps(k, 2), 'wo', ...
                'MarkerSize', 8 + 12*comps(k,3)/max(comps(:,3)), ...
                'LineWidth', 1.5);
        end
    end
    xlabel('Range (m)');
    ylabel('Cross-Range (m)');
    title(sprintf('PSF+CLEAN Δθ=%.1f° (%d centers)', ...
        aperture_info{idx}.delta_theta_deg, size(comps,1)));
    caxis(caxis_range);
    colorbar;
    hold off;
end

sgtitle(sprintf('Method 2: PSF Superposition SAR Imaging | Model: %s | θ_0=%.0f°', ...
    modelName, theta_center_deg), 'FontSize', 13, 'FontWeight', 'bold');

%% --- 图3: FFT vs PSF 并排对比 ---
figure(3);
clf;
set(gcf, 'Name', 'FFT vs PSF Comparison (FFT vs PSF对比)', ...
         'NumberTitle', 'off', 'Position', [150, 80, 350*N_valid, 600]);

for i = 1:N_valid
    idx = valid_idx(i);
    r_ax = results_fft{idx}.range_axis;
    x_ax = results_fft{idx}.crossrange_axis;

    % FFT
    subplot(3, N_valid, i);
    imagesc(r_ax, x_ax, results_fft{idx}.dB);
    colormap('jet');
    axis xy;
    xlabel('Range (m)'); ylabel('Cross-Range (m)');
    title(sprintf('FFT Δθ=%.1f°', aperture_info{idx}.delta_theta_deg));
    caxis(caxis_range); colorbar;

    % PSF
    subplot(3, N_valid, N_valid + i);
    imagesc(r_ax, x_ax, results_psf{idx}.dB);
    colormap('jet');
    axis xy;
    xlabel('Range (m)'); ylabel('Cross-Range (m)');
    title(sprintf('PSF Δθ=%.1f°', aperture_info{idx}.delta_theta_deg));
    caxis(caxis_range); colorbar;

    % Difference
    diff_dB = results_fft{idx}.dB - results_psf{idx}.dB;
    subplot(3, N_valid, 2*N_valid + i);
    imagesc(r_ax, x_ax, diff_dB);
    colormap('jet');
    axis xy;
    xlabel('Range (m)'); ylabel('Cross-Range (m)');
    title(sprintf('Diff (FFT-PSF) Δθ=%.1f°', aperture_info{idx}.delta_theta_deg));
    caxis([-10, 10]); colorbar;
end

sgtitle(sprintf('FFT vs PSF Comparison | Model: %s | θ_0=%.0f°', ...
    modelName, theta_center_deg), 'FontSize', 13, 'FontWeight', 'bold');

%% --- 图4: 指标对比曲线 ---
if N_valid >= 1
    figure(4);
    clf;
    set(gcf, 'Name', 'Metrics Comparison (指标对比)', ...
             'NumberTitle', 'off', 'Position', [200, 100, 800, 600]);

    delta_theta_vals = [metrics.delta_theta];

    % 相关系数对比
    subplot(2, 2, 1);
    plot(delta_theta_vals, [metrics.cc], 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Angular Aperture Δθ (deg)');
    ylabel('Correlation Coefficient (CC)');
    title('CC: FFT vs PSF');
    grid on;
    xlim([min(delta_theta_vals)*0.8, max(delta_theta_vals)*1.2]);
    ylim([0, 1]);

    % 平均幅度对比 (FFT vs PSF)
    subplot(2, 2, 2);
    plot(delta_theta_vals, [metrics.amp_fft], 'bs-', 'LineWidth', 1.5, 'MarkerSize', 8);
    hold on;
    plot(delta_theta_vals, [metrics.amp_psf], 'r^-', 'LineWidth', 1.5, 'MarkerSize', 8);
    hold off;
    xlabel('Angular Aperture Δθ (deg)');
    ylabel('Mean Amplitude');
    title('Mean Amplitude: FFT vs PSF');
    legend('FFT', 'PSF (CLEAN)', 'Location', 'best');
    grid on;
    xlim([min(delta_theta_vals)*0.8, max(delta_theta_vals)*1.2]);

    % 展宽宽度对比
    subplot(2, 2, 3);
    plot(delta_theta_vals, [metrics.bw_fft], 'bs-', 'LineWidth', 1.5, 'MarkerSize', 8);
    hold on;
    plot(delta_theta_vals, [metrics.bw_psf], 'r^-', 'LineWidth', 1.5, 'MarkerSize', 8);
    hold off;
    xlabel('Angular Aperture Δθ (deg)');
    ylabel('Beamwidth (bins)');
    title('Azimuth Beamwidth (-3dB)');
    legend('FFT', 'PSF', 'Location', 'best');
    grid on;

    % 方位分辨率对比
    subplot(2, 2, 4);
    delta_x_vals = [metrics.delta_x];
    plot(delta_theta_vals, delta_x_vals, 'go-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Angular Aperture Δθ (deg)');
    ylabel('Δx Cross-Range Resolution (m)');
    title('Theoretical Cross-Range Resolution');
    grid on;
    xlim([min(delta_theta_vals)*0.8, max(delta_theta_vals)*1.2]);

    sgtitle(sprintf('Metrics vs Angular Aperture | Model: %s | θ_0=%.0f°', ...
        modelName, theta_center_deg), 'FontSize', 13, 'FontWeight', 'bold');
end

%% --- 图5: CLEAN 散射中心分布（最大孔径）---
if N_valid >= 1 && ~isempty(results_clean{valid_idx(end)})
    figure(5);
    clf;
    set(gcf, 'Name', 'CLEAN Scattering Centers (散射中心分布)', ...
             'NumberTitle', 'off', 'Position', [250, 120, 700, 550]);

    % 使用最大孔径
    idx_last = valid_idx(end);
    r_ax = results_fft{idx_last}.range_axis;
    x_ax = results_fft{idx_last}.crossrange_axis;
    comps = results_clean{idx_last};

    % 背景：FFT 图像
    imagesc(r_ax, x_ax, results_fft{idx_last}.dB);
    colormap('jet');
    axis xy;
    hold on;

    % 标注 CLEAN 散射中心
    if ~isempty(comps)
        amplitudes = comps(:, 3);
        amp_norm = amplitudes / max(amplitudes);
        for k = 1:size(comps, 1)
            sz = 20 + 80 * amp_norm(k);
            scatter(comps(k, 1), comps(k, 2), sz, 'wo', 'LineWidth', 1.5);
            text(comps(k, 1) + 0.02, comps(k, 2) + 0.02, num2str(k), ...
                'Color', 'w', 'FontSize', 8, 'FontWeight', 'bold');
        end
    end
    hold off;

    xlabel('Range (m)');
    ylabel('Cross-Range (m)');
    title(sprintf('CLEAN Scattering Centers (%d) Overlaid on FFT Image | Δθ=%.1f°', ...
        size(comps, 1), aperture_info{idx_last}.delta_theta_deg), ...
        'FontSize', 13, 'FontWeight', 'bold');
    caxis(caxis_range);
    colorbar;
end

%% ========================================================================
%% 6. 保存结果
%% ========================================================================
fprintf('\nSaving results...\n');

% Create timestamped result directory
[resultDir, nowStr] = createResultDir('main_sar_imaging');

% 保存各图
for figNum = 1:5
    if ishandle(figNum)
        figFile = fullfile(resultDir, sprintf('sar_fig%d_%s.png', figNum, nowStr));
        saveas(figNum, figFile);
        fprintf('  Figure %d: %s\n', figNum, figFile);
    end
end

% 保存 MAT 数据
sarDataFile = fullfile(resultDir, ['sar_imaging_' nowStr '.mat']);
save(sarDataFile, 'results_fft', 'results_psf', 'results_clean', ...
    'aperture_info', 'metrics', 'valid_idx', 'theta_center_deg', ...
    'apertures_deg', 'f_center', 'lambda_c', 'B', 'delta_r', 'c', ...
    'windowType', 'cleanThreshold', 'cleanMaxIter', '-v7.3');
fprintf('  MAT data: %s\n', sarDataFile);

% 导出文本结果
txtFile = fullfile(resultDir, ['sar_metrics_' nowStr '.txt']);
fid = fopen(txtFile, 'w');
fprintf(fid, '# SAR Imaging Metrics Summary\n');
fprintf(fid, '# Generated: %s\n', nowStr);
fprintf(fid, '# Model: %s, Center angle: %.1f°\n', modelName, theta_center_deg);
fprintf(fid, '# f_c: %.2f GHz, λ_c: %.4f m, B: %.2f GHz, Δr: %.4f m\n', ...
    f_center/1e9, lambda_c, B/1e9, delta_r);
fprintf(fid, '#\n');
fprintf(fid, '# Δθ(°)  Δx(m)  N_θ  CC      SSIM    BW_FFT  BW_PSF\n');
fprintf(fid, '# ------ ------ ---- ------- ------- ------- -------\n');
for i = 1:N_valid
    fprintf(fid, '  %5.1f  %5.2f  %4d  %7.4f %7.4f %7.2f %7.2f\n', ...
        metrics(i).delta_theta, metrics(i).delta_x, ...
        aperture_info{valid_idx(i)}.N_theta_sub, ...
        metrics(i).cc, metrics(i).ssim, ...
        metrics(i).bw_fft, metrics(i).bw_psf);
end
fclose(fid);
fprintf('  Text summary: %s\n', txtFile);

%% ========================================================================
%% 7. 显示摘要
%% ========================================================================
fprintf('\n========================================\n');
fprintf('  SAR Imaging Complete\n');
fprintf('========================================\n');
fprintf('  Model:               %s\n', modelName);
fprintf('  Center frequency:    %.2f GHz\n', f_center/1e9);
fprintf('  Wavelength:          %.4f m\n', lambda_c);
fprintf('  Bandwidth:           %.2f GHz\n', B/1e9);
fprintf('  Range resolution:    %.4f m\n', delta_r);
fprintf('  Center angle θ0:     %.1f°\n', theta_center_deg);
fprintf('  Sub-apertures:       %s\n', strjoin(aperture_labels, ', '));
fprintf('  Window:              %s\n', windowType);
fprintf('\n  Aperture Comparison:\n');
fprintf('  %-8s %-8s %-8s %-8s %-12s %-12s\n', ...
    'Δθ(°)', 'Δx(m)', 'CC', 'SSIM', 'BW_FFT(bin)', 'BW_PSF(bin)');
fprintf('  %-8s %-8s %-8s %-8s %-12s %-12s\n', ...
    '------', '------', '------', '------', '-----------', '-----------');
for i = 1:N_valid
    fprintf('  %-8.1f %-8.2f %-8.4f %-8.4f %-12.2f %-12.2f\n', ...
        metrics(i).delta_theta, metrics(i).delta_x, ...
        metrics(i).cc, metrics(i).ssim, ...
        metrics(i).bw_fft, metrics(i).bw_psf);
end
fprintf('========================================\n');

end  % main_sar_imaging

%% ========================================================================
%% 辅助函数
%% ========================================================================

% -------------------------------------------------------------------------
function [S_sub, theta_sub, N_theta_sub] = extractSubAperture(S_complex, ...
    theta_vec, theta_min, theta_max)
% 从全数据中提取指定角度范围的子孔径数据
%
% Input:
%   S_complex  - 复散射场矩阵 (N_angles × N_f)
%   theta_vec  - 角度向量 (N_angles × 1)
%   theta_min  - 角度下限 (deg)
%   theta_max  - 角度上限 (deg)
% Output:
%   S_sub      - 子孔径复散射场 (N_theta_sub × N_f)
%   theta_sub  - 角度向量 (N_theta_sub × 1)
%   N_theta_sub - 子孔径角度数

    % 找到在 [theta_min, theta_max] 内的角度索引
    angle_indices = find(theta_vec >= theta_min & theta_vec <= theta_max);

    if isempty(angle_indices)
        % 宽松匹配：找最近的 N 个角度
        [~, sort_idx] = sort(abs(theta_vec - (theta_min + theta_max)/2));
        half_N = max(1, floor(length(sort_idx) * (theta_max - theta_min) / ...
                                (theta_vec(end) - theta_vec(1)) / 4));
        angle_indices = sort_idx(1:min(2*half_N, length(sort_idx)));
        angle_indices = sort(angle_indices);
        warning('  No exact angle match, using %d nearest angles.', length(angle_indices));
    end

    S_sub = S_complex(angle_indices, :);
    theta_sub = theta_vec(angle_indices);
    N_theta_sub = length(theta_sub);
end

% -------------------------------------------------------------------------
function [W_2d, w_f, w_theta] = create2DWindow(N_f, N_theta, windowType)
% 创建二维可分离窗函数
%
% Input:
%   N_f        - 频率点数
%   N_theta    - 角度点数
%   windowType - 窗类型: 'hamming', 'hanning', 'blackman', 'kaiser', 'none'
% Output:
%   W_2d       - 二维窗矩阵 (N_theta × N_f)
%   w_f        - 频率维窗向量 (1 × N_f)
%   w_theta    - 角度维窗向量 (N_theta × 1)

    switch lower(windowType)
        case 'hamming'
            w_f = hamming(N_f)';
            w_theta = hamming(N_theta);
        case 'hanning'
            w_f = hanning(N_f)';
            w_theta = hanning(N_theta);
        case 'blackman'
            w_f = blackman(N_f)';
            w_theta = blackman(N_theta);
        case 'kaiser'
            w_f = kaiser(N_f, 2.5)';
            w_theta = kaiser(N_theta, 2.5);
        case 'none'
            w_f = ones(1, N_f);
            w_theta = ones(N_theta, 1);
        otherwise
            warning('Unknown window type "%s", using Hamming.', windowType);
            w_f = hamming(N_f)';
            w_theta = hamming(N_theta);
    end

    % 确保正确的维度
    w_f = w_f(:)';
    w_theta = w_theta(:);

    % 二维窗: 外积 w_theta * w_f
    W_2d = w_theta * w_f;
end

% -------------------------------------------------------------------------
function [components, residual] = cleanAlgorithm(I_complex, range_axis, ...
    crossrange_axis, psf_func, threshold, maxIter, gain)
% CLEAN 迭代算法：从复图像中提取点散射中心
%
% Input:
%   I_complex       - 复图像 (N_x_pad × N_r_pad)
%   range_axis      - 距离轴 (m)
%   crossrange_axis - 方位轴 (m)
%   psf_func        - PSF 函数句柄: h(r, x)
%   threshold       - 能量阈值 (相对于初始峰值)
%   maxIter         - 最大迭代次数
%   gain            - CLEAN 增益因子 (0 < gain <= 1)
% Output:
%   components      - K × 4 矩阵: [r_k, x_k, A_k, phi_k] (位置+幅度+相位)
%   residual        - 残差图像

    I_mag = abs(I_complex);
    N_r = length(range_axis);
    N_x = length(crossrange_axis);

    % 初始化
    residual = I_complex;
    residual_mag = I_mag;
    components = [];

    % 初始峰值作为阈值参考
    initial_peak = max(residual_mag(:));
    stop_threshold = threshold * initial_peak;

    fprintf('    CLEAN: initial peak=%.4e, stop threshold=%.4e\n', initial_peak, stop_threshold);

    for iter = 1:maxIter
        % 寻找当前残差的峰值
        [peak_val, peak_idx] = max(residual_mag(:));

        if peak_val < stop_threshold
            break;
        end

        % 峰值位置
        [peak_row, peak_col] = ind2sub([N_x, N_r], peak_idx);
        r_peak = range_axis(peak_col);
        x_peak = crossrange_axis(peak_row);

        % 峰值复值
        complex_peak = residual(peak_row, peak_col);
        amplitude = abs(complex_peak);
        phase_val = angle(complex_peak);

        % 记录散射中心
        components = [components; r_peak, x_peak, amplitude, phase_val]; %#ok<AGROW>

        % 构建当前点处的 PSF（直接求值，移位到峰值位置）
        psf_shifted = buildShiftedPSF(range_axis, crossrange_axis, ...
            psf_func, r_peak, x_peak);

        % 从残差中减去 (gain * complex_peak * PSF)
        residual = residual - gain * complex_peak * psf_shifted;
        residual_mag = abs(residual);
    end

    fprintf('    CLEAN converged after %d iterations, %.2f%% residual energy\n', ...
        iter, 100 * max(residual_mag(:)) / initial_peak);
end

% -------------------------------------------------------------------------
function [I_psf, I_psf_mag] = buildPSFImage(components, range_axis, ...
    crossrange_axis, psf_func)
% 基于 CLEAN 提取的散射中心，通过 PSF 叠加构建图像
%
% Input:
%   components     - K × 4 矩阵: [r_k, x_k, A_k, phi_k]
%   range_axis     - 距离轴 (m)
%   crossrange_axis - 方位轴 (m)
%   psf_func       - PSF 函数句柄
% Output:
%   I_psf          - 复图像
%   I_psf_mag      - 幅度图像

    N_r = length(range_axis);
    N_x = length(crossrange_axis);

    I_psf = zeros(N_x, N_r);

    if isempty(components)
        I_psf_mag = abs(I_psf);
        return;
    end

    % 对每个散射中心进行 PSF 叠加（直接求值法）
    for k = 1:size(components, 1)
        r_k = components(k, 1);
        x_k = components(k, 2);
        A_k = components(k, 3);
        phi_k = components(k, 4);

        % 复增益
        gain_c = A_k * exp(1j * phi_k);

        % 移位 PSF
        psf_shifted = buildShiftedPSF(range_axis, crossrange_axis, ...
            psf_func, r_k, x_k);

        % 叠加
        I_psf = I_psf + gain_c * psf_shifted;
    end

    I_psf_mag = abs(I_psf);
end

% -------------------------------------------------------------------------
function psf_shifted = buildShiftedPSF(range_axis, crossrange_axis, ...
    psf_func, r_center, x_center)
% 构建移位到指定散射中心位置的 PSF（矢量化版本）
%
% Input:
%   range_axis      - 距离轴 (m)
%   crossrange_axis - 方位轴 (m)
%   psf_func        - PSF 函数句柄: h(r, x)
%   r_center        - 散射中心距离位置 (m)
%   x_center        - 散射中心方位位置 (m)
% Output:
%   psf_shifted     - 移位后的 PSF (N_x × N_r)

    [R_mesh, X_mesh] = meshgrid(range_axis, crossrange_axis);
    psf_shifted = psf_func(R_mesh - r_center, X_mesh - x_center);
end

% -------------------------------------------------------------------------
function cc = computeCorrelationCoefficient(X, Y)
% 计算两幅图像的相关系数 (Correlation Coefficient)
%
%   ρ = Σ((X-μ_X)·(Y-μ_Y)) / √(Σ(X-μ_X)² · Σ(Y-μ_Y)²)

    X = double(X(:));
    Y = double(Y(:));

    X_mean = mean(X);
    Y_mean = mean(Y);

    numerator = sum((X - X_mean) .* (Y - Y_mean));
    denominator = sqrt(sum((X - X_mean).^2) * sum((Y - Y_mean).^2));

    if denominator < 1e-15
        cc = 0;
    else
        cc = numerator / denominator;
    end
end

% -------------------------------------------------------------------------
function ssim_val = computeSSIM(X, Y)
% 计算结构相似度 (SSIM)
%
%   SSIM(x,y) = (2μ_x μ_y + C1)(2σ_xy + C2) /
%               ((μ_x² + μ_y² + C1)(σ_x² + σ_y² + C2))
%
% 使用高斯滑动窗口的均值；小图像回退到全局SSIM

    X = double(X);
    Y = double(Y);

    % 动态范围
    L_dr = max(max(X(:)), max(Y(:))) - min(min(X(:)), min(Y(:)));
    if L_dr < 1e-10, L_dr = 1; end

    % 常数
    K1 = 0.01;
    K2 = 0.03;
    C1 = (K1 * L_dr)^2;
    C2 = (K2 * L_dr)^2;

    min_dim = min(size(X));

    if min_dim >= 8
        % 使用局部滑动窗口 SSIM
        win_size = min(8, min_dim);
        sigma_win = 1.5;
        [gx, gy] = meshgrid(-(win_size-1)/2:(win_size-1)/2);
        gauss_win = exp(-(gx.^2 + gy.^2) / (2 * sigma_win^2));
        gauss_win = gauss_win / sum(gauss_win(:));

        mu_x = conv2(X, gauss_win, 'valid');
        mu_y = conv2(Y, gauss_win, 'valid');

        if isempty(mu_x)
            % 回退到全局
            mu_x_global = mean(X(:));
            mu_y_global = mean(Y(:));
            sigma_x_sq = mean((X(:) - mu_x_global).^2);
            sigma_y_sq = mean((Y(:) - mu_y_global).^2);
            sigma_xy = mean((X(:) - mu_x_global) .* (Y(:) - mu_y_global));

            numerator = (2 * mu_x_global * mu_y_global + C1) * (2 * sigma_xy + C2);
            denominator = (mu_x_global^2 + mu_y_global^2 + C1) * (sigma_x_sq + sigma_y_sq + C2);
            ssim_val = numerator / denominator;
        else
            sigma_x_sq = conv2(X.^2, gauss_win, 'valid') - mu_x.^2;
            sigma_y_sq = conv2(Y.^2, gauss_win, 'valid') - mu_y.^2;
            sigma_xy = conv2(X .* Y, gauss_win, 'valid') - mu_x .* mu_y;

            numerator = (2 * mu_x .* mu_y + C1) .* (2 * sigma_xy + C2);
            denominator = (mu_x.^2 + mu_y.^2 + C1) .* (sigma_x_sq + sigma_y_sq + C2);

            ssim_map = numerator ./ denominator;
            ssim_val = mean(ssim_map(:));
        end
    else
        % 小图像：使用全局 SSIM
        mu_x_global = mean(X(:));
        mu_y_global = mean(Y(:));
        sigma_x_sq = mean((X(:) - mu_x_global).^2);
        sigma_y_sq = mean((Y(:) - mu_y_global).^2);
        sigma_xy = mean((X(:) - mu_x_global) .* (Y(:) - mu_y_global));

        numerator = (2 * mu_x_global * mu_y_global + C1) * (2 * sigma_xy + C2);
        denominator = (mu_x_global^2 + mu_y_global^2 + C1) * (sigma_x_sq + sigma_y_sq + C2);
        ssim_val = numerator / denominator;
    end

    if isnan(ssim_val) || ssim_val < 0
        ssim_val = 0;
    elseif ssim_val > 1
        ssim_val = 1;
    end
end

% -------------------------------------------------------------------------
function bw = measureBeamwidth(I_dB)
% 测量方位向 -3dB 主瓣宽度
%
% Input:
%   I_dB  - dB 归一化的图像 (N_x × N_r)
% Output:
%   bw    - 平均 -3dB 宽度（像素数）

    % 对每列（每个距离单元）沿方位向找 -3dB 宽度
    N_x = size(I_dB, 1);
    N_r = size(I_dB, 2);

    widths = [];

    for col = 1:N_r
        profile = I_dB(:, col);
        peak_val = max(profile);

        if peak_val < -30  % 跳过没有明显峰值的列
            continue;
        end

        % 找峰值位置
        [~, peak_idx] = max(profile);

        % 从峰值向两侧找 -3dB 点
        thresh = peak_val - 3;

        % 左侧
        left_idx = peak_idx;
        while left_idx > 1 && profile(left_idx) > thresh
            left_idx = left_idx - 1;
        end

        % 右侧
        right_idx = peak_idx;
        while right_idx < N_x && profile(right_idx) > thresh
            right_idx = right_idx + 1;
        end

        width = right_idx - left_idx;
        if width > 0 && width < N_x
            widths = [widths, width]; %#ok<AGROW>
        end
    end

    if isempty(widths)
        bw = 0;
    else
        bw = mean(widths);
    end
end

% -------------------------------------------------------------------------
function y = sincf(x)
% Sinc 函数: sin(pi*x) / (pi*x)
% MATLAB 没有内置 sinc 函数，此处定义。
% 在 x=0 处返回 1（极限值）。
%
% Usage:
%   y = sincf(x)  返回 sin(pi*x)/(pi*x)，x 可以是标量、向量或矩阵

    % 处理 x=0 的奇点
    y = ones(size(x));
    nonzero = abs(x) > 1e-15;
    y(nonzero) = sin(pi * x(nonzero)) ./ (pi * x(nonzero));
end
