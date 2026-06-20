function main_scattering_center_modeling(varargin)
% MAIN_SCATTERING_CENTER_MODELING  基于PO结果的散射中心建模
%
%   从 PO (Physical Optics) 宽带复散射场数据出发，自动提取并建立目标的
%   散射中心参数化模型。完整实现了"基于PO结果的散射中心建模方法.md"中
%   描述的六步流程：
%
%     步骤1: 生成 HRRP 历程图（IFFT沿频率维）
%     步骤2: 峰值检测与候选提取（CLEAN算法）
%     步骤3: 参数初步估计（幅度A、频率依赖α、角度依赖特性）
%     步骤4: 散射中心分类（局部型/分布型/滑动型）
%     步骤5: 参数优化（非线性最小二乘）
%     步骤6: 模型验证（合成回波 vs PO回波对比）
%
%   GTD (Geometrical Theory of Diffraction) 散射中心模型:
%
%     局部型:   S_k(f,φ) = A_k * (j*f/f_c)^α_k * exp(-j*4πf/c * R_k)
%     分布型:   S_k(f,φ) = A_k * (j*f/f_c)^α_k * exp(-j*4πf/c * R_k)
%                          * sinc(2πf/c * L_k * sin(φ-φ̄_k))
%     滑动型:   S_k(f,φ) = A_k * (j*f/f_c)^α_k * exp(-j*4πf/c * R_k(φ))
%                          R_k(φ) = R_k0 + dR/dφ * (φ-φ̄_k)
%
%   输入文件:
%     results/wideband_scattering_<timestamp>.mat  (由 main_wideband_scattering 生成)
%
%   Usage:
%     >> main_scattering_center_modeling                % 自动查找最新数据
%     >> main_scattering_center_modeling('file.mat')    % 指定数据文件
%     >> main_scattering_center_modeling(..., 'MaxCenters', 15, 'CleanGain', 0.7)

%% ========================================================================
%% 0. 参数解析
%% ========================================================================
addpath('lib');

p = inputParser;
p.addOptional('DataFile', '', @(x) ischar(x) || isstring(x) || isempty(x));
p.addParameter('MaxCenters', 15, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('CleanGain', 0.7, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
p.addParameter('CleanThreshold', 0.02, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('ZeroPadFactor', 4, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('WindowType', 'hamming', @(x) ischar(x) || isstring(x));
p.addParameter('SlidingThreshold', 1.0, @(x) isnumeric(x) && isscalar(x));
p.addParameter('DistributedWidthFactor', 1.5, @(x) isnumeric(x) && isscalar(x));
p.addParameter('OptimizeParams', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('ShowDetailed', true, @(x) islogical(x) || isnumeric(x));

if ~isempty(varargin) && (ischar(varargin{1}) || isstring(varargin{1})) ...
        && any(strcmp(char(varargin{1}), {'MaxCenters','CleanGain','CleanThreshold',...
        'ZeroPadFactor','WindowType','SlidingThreshold','DistributedWidthFactor',...
        'OptimizeParams','ShowDetailed'}))
    varargin = [{''}, varargin];
end
p.parse(varargin{:});

dataFile              = p.Results.DataFile;
K_max                 = p.Results.MaxCenters;
cleanGain             = p.Results.CleanGain;
cleanThreshold        = p.Results.CleanThreshold;
zeroPadFactor         = p.Results.ZeroPadFactor;
windowType            = p.Results.WindowType;
slidingThreshold      = p.Results.SlidingThreshold;
distributedWidthFactor = p.Results.DistributedWidthFactor;
doOptimize            = p.Results.OptimizeParams;
showDetailed          = p.Results.ShowDetailed;

fprintf('========================================\n');
fprintf('  Scattering Center Modeling from PO\n');
fprintf('  (基于PO结果的散射中心建模)\n');
fprintf('========================================\n\n');
fprintf('Parameters:\n');
fprintf('  Max centers:    %d\n', K_max);
fprintf('  CLEAN gain:     %.2f\n', cleanGain);
fprintf('  CLEAN threshold: %.3f\n', cleanThreshold);
fprintf('  Zero-pad factor: %dx\n', zeroPadFactor);
fprintf('  Window:         %s\n', windowType);
fprintf('  Optimization:   %s\n', ternary(doOptimize, 'on', 'off'));

%% ========================================================================
%% 步骤0: 加载PO宽带散射数据
%% ========================================================================
fprintf('\n========================================\n');
fprintf('  Loading PO Wideband Scattering Data\n');
fprintf('========================================\n');

if isempty(dataFile)
    dataFile = findLatestResultFile('wideband_scattering_*.mat');
    if isempty(dataFile)
        error('No wideband_scattering_*.mat found in results/.\nRun main_wideband_scattering first.');
    end
end

fprintf('  File: %s\n', dataFile);
data = load(dataFile);

% 验证必要字段
requiredFields = {'S_complex', 'freq_array', 'theta_array', 'phi_array', ...
                  'c', 'B', 'delta_r', 'ip', 'it', 'N_f', 'N_angles_total'};
for k = 1:length(requiredFields)
    if ~isfield(data, requiredFields{k})
        error('Missing required field: %s', requiredFields{k});
    end
end

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

% 角度步长与范围
if isfield(data, 'delt'), delt = data.delt; else delt = 1; end
if isfield(data, 'delp'), delp = data.delp; else delp = 1; end
if isfield(data, 'tstart'), tstart = data.tstart; else tstart = 0; end
if isfield(data, 'tstop'),  tstop  = data.tstop;  else tstop  = 180; end
if isfield(data, 'pstart'), pstart = data.pstart; else pstart = 0; end

% 中心频率（用于GTD模型参考频率）
f_c = (freq_array(1) + freq_array(end)) / 2;
lambda_c = c / f_c;
freq_step = freq_array(2) - freq_array(1);

% 构建角度向量（与main_range_profile.m一致）
if ip > 1 && it == 1
    angle_vec = phi_array(:, 1);
    angleLabel = '\phi';
    angleUnit = 'deg';
elseif it > 1 && ip == 1
    angle_vec = theta_array(1, :)';
    angleLabel = '\theta';
    angleUnit = 'deg';
else
    angle_vec = (1:N_angles_total)';
    angleLabel = 'idx';
    angleUnit = '';
end
N_angles = length(angle_vec);

fprintf('  Model:          %s\n', modelName);
fprintf('  Frequency:      %.2f - %.2f GHz (%d points)\n', freq_array(1)/1e9, freq_array(end)/1e9, N_f);
fprintf('  Bandwidth B:    %.2f GHz\n', B/1e9);
fprintf('  Center freq fc: %.2f GHz\n', f_c/1e9);
fprintf('  Range res Δr:   %.4f m\n', delta_r);
fprintf('  Angles:         %d points, %.1f° to %.1f°\n', N_angles, angle_vec(1), angle_vec(end));

%% ========================================================================
%% 步骤1: 生成 HRRP 历程图
%% ========================================================================
fprintf('\n========================================\n');
fprintf('  Step 1: Generating HRRP History Map\n');
fprintf('  (一维距离像历程图)\n');
fprintf('========================================\n');

% 配置窗函数
switch lower(windowType)
    case 'hamming'
        W_freq = hamming(N_f)';
    case 'hanning'
        W_freq = hanning(N_f)';
    case 'blackman'
        W_freq = blackman(N_f)';
    case 'chebwin'
        W_freq = chebwin(N_f, 60)';
    case 'kaiser'
        W_freq = kaiser(N_f, 2.5)';
    case 'none'
        W_freq = ones(1, N_f);
    otherwise
        warning('Unknown window type "%s", using Hamming.', windowType);
        W_freq = hamming(N_f)';
end
W_freq = W_freq(:)';

% 补零FFT
N_fft = N_f * zeroPadFactor;
if mod(N_fft, 2) ~= 0
    N_fft = N_fft + 1;
end

% 距离轴标定
delta_r_fft = delta_r * (N_f / N_fft);
range_axis = (-floor(N_fft/2) : floor(N_fft/2)-1) * delta_r_fft;

fprintf('  FFT points: %d (zero-pad: %dx)\n', N_fft, zeroPadFactor);
fprintf('  Range bin spacing: %.4f m\n', delta_r_fft);
fprintf('  Range axis: %.2f ~ %.2f m\n', range_axis(1), range_axis(end));

% 预计算所有角度的 HRRP
fprintf('  Computing HRRP for all %d angles...\n', N_angles);
tic;

% HRRP 复数 (N_angles × N_fft) — 保留复数用于相位分析
HRRP_complex = zeros(N_angles, N_fft);

for i_ang = 1:N_angles
    Es_f = S_complex(i_ang, :);  % 1 × N_f
    Es_windowed = W_freq .* Es_f;
    Es_padded = [Es_windowed, zeros(1, N_fft - N_f)];
    h_complex = ifft(Es_padded, N_fft);
    h_complex = ifftshift(h_complex);
    HRRP_complex(i_ang, :) = h_complex;
end

% HRRP 幅度 (dB) — 用于峰值检测
HRRP_mag = abs(HRRP_complex);
epsilon_dB = 1e-10;
HRRP_dB = 20 * log10(HRRP_mag + epsilon_dB);

elapsed = toc;
fprintf('  HRRP computation: %.2f s (%.3f s/angle)\n', elapsed, elapsed/N_angles);

%% ========================================================================
%% 步骤2-3: CLEAN 算法 — 峰值检测与参数估计
%% ========================================================================
fprintf('\n========================================\n');
fprintf('  Step 2-3: CLEAN Peak Detection & Parameter Estimation\n');
fprintf('  (CLEAN峰值检测与参数估计)\n');
fprintf('========================================\n');

% 初始化
residual = S_complex;           % 频域残差 (N_angles × N_f)
centers = struct([]);           % 散射中心参数列表
K = 0;                         % 已提取的散射中心数

% 初始信号能量（用于CLEAN终止判据）
initial_energy = sum(abs(S_complex(:)).^2);

% 用于存储每次迭代的信息
clean_history = [];

fprintf('  Initial signal energy: %.4e\n', initial_energy);
fprintf('  CLEAN gain: %.2f, threshold: %.3f, max centers: %d\n', ...
    cleanGain, cleanThreshold, K_max);

for iter = 1:K_max
    %% --- 2a. 对当前残差求HRRP ---
    residual_HRRP = zeros(N_angles, N_fft);
    for i_ang = 1:N_angles
        Es_f = residual(i_ang, :);
        Es_windowed = W_freq .* Es_f;
        Es_padded = [Es_windowed, zeros(1, N_fft - N_f)];
        h = ifft(Es_padded, N_fft);
        h = ifftshift(h);
        residual_HRRP(i_ang, :) = h;
    end
    residual_HRRP_mag = abs(residual_HRRP);

    %% --- 2b. 找到最强峰值 ---
    [peak_val, peak_lin_idx] = max(residual_HRRP_mag(:));
    [peak_ang_idx, peak_r_idx] = ind2sub([N_angles, N_fft], peak_lin_idx);

    R_k = range_axis(peak_r_idx);
    phi_k = angle_vec(peak_ang_idx);

    residual_energy = sum(abs(residual(:)).^2);
    energy_ratio = residual_energy / initial_energy;

    if energy_ratio < cleanThreshold
        fprintf('  CLEAN converged: residual energy %.2e (%.2f%% of initial)\n', ...
            residual_energy, 100*energy_ratio);
        break;
    end

    fprintf('\n  --- Iteration %d ---\n', iter);
    fprintf('  Peak: R=%.4f m, %s=%.1f°, |h|=%.4e (%.1f dB)\n', ...
        R_k, angleLabel, phi_k, peak_val, 20*log10(peak_val));

    %% --- 3a. 估计频率依赖指数 α_k ---
    % 取峰值角度处的频域数据，相位解调后做 log|S| vs log(f) 线性拟合
    S_peak_f = residual(peak_ang_idx, :);  % 1 × N_f

    % 相位解调：去除距离相位 exp(-j*4πf/c * R_k)
    phase_demod = exp(1j * (4*pi/c) * R_k * freq_array);
    S_demod = S_peak_f .* phase_demod;

    % 取幅度
    S_demod_mag = abs(S_demod);

    % 有效频率范围（排除幅度过低的频点）
    valid_idx = S_demod_mag > max(S_demod_mag) * 1e-4;
    if sum(valid_idx) < 5
        valid_idx = true(size(S_demod_mag));
    end

    log_f = log(freq_array(valid_idx) / f_c);
    log_S = log(S_demod_mag(valid_idx));

    % 线性拟合: log|S| = α * log(f/f_c) + const
    p = polyfit(log_f, log_S, 1);
    alpha_k = p(1);

    % 拟合质量
    log_S_fit = polyval(p, log_f);
    R2 = 1 - sum((log_S - log_S_fit).^2) / sum((log_S - mean(log_S)).^2);

    fprintf('  Estimated α = %.3f (R² = %.3f)\n', alpha_k, R2);

    %% --- 3b. 估计幅度 A_k ---
    A_k = peak_val;

    %% --- 3c. 跟踪峰值位置随角度的变化 ---
    % 在峰值角度附近跟踪R_peak，判断是否为滑动型
    angle_half_span = max(5, floor(N_angles * 0.1));  % 搜索范围：±10%的角度
    ang_start = max(1, peak_ang_idx - angle_half_span);
    ang_end = min(N_angles, peak_ang_idx + angle_half_span);
    N_track = ang_end - ang_start + 1;

    R_track = zeros(N_track, 1);
    amp_track = zeros(N_track, 1);
    phi_track = angle_vec(ang_start:ang_end);

    for i_t = 1:N_track
        i_ang = ang_start + i_t - 1;
        [local_peak, local_r_idx] = max(residual_HRRP_mag(i_ang, :));
        R_track(i_t) = range_axis(local_r_idx);
        amp_track(i_t) = local_peak;
    end

    % 估计位置变化率 dR/dφ
    if N_track >= 3
        phi_rad = deg2rad(phi_track);
        p_R = polyfit(phi_rad, R_track, 1);
        dR_dphi = p_R(1);  % m/rad

        % R² of linear fit
        R_fit = polyval(p_R, phi_rad);
        R2_sliding = 1 - sum((R_track - R_fit).^2) / sum((R_track - mean(R_track)).^2);
    else
        dR_dphi = 0;
        R2_sliding = 0;
    end

    % 位置变化量（在一个距离分辨单元中）
    delta_R_range = max(R_track) - min(R_track);
    delta_R_norm = delta_R_range / delta_r;  % 归一化到分辨率单元

    fprintf('  R tracking: ΔR=%.4f m (%.1f×Δr), dR/dφ=%.3f m/rad (R²=%.3f)\n', ...
        delta_R_range, delta_R_norm, dR_dphi, R2_sliding);

    %% --- 3d. 估计角度展宽和有效长度 L_k ---
    % 对角度幅度剖面做 sinc 拟合估计 L_k
    amp_norm = amp_track / max(amp_track);

    % 找-6 dB角度宽度
    above_half = amp_norm > 0.5;
    if sum(above_half) >= 2
        half_indices = find(above_half);
        angle_width_deg = phi_track(half_indices(end)) - phi_track(half_indices(1));
    else
        angle_width_deg = 0;
    end

    % 从sinc宽度估计L_k: sinc主瓣零点宽度 = c/(f*L*cos(φ-φ̄))
    % -6dB宽度 ≈ 0.6 * 主瓣宽度
    if angle_width_deg > 0 && angle_width_deg < 180
        angle_width_rad = deg2rad(angle_width_deg);
        % sinc(2πf/c * L * sin(φ-φ̄)) 的-6dB点
        % 近似: L ≈ 0.6 * c / (f_c * angle_width_rad) (小角度)
        L_k_est = 0.6 * c / (f_c * angle_width_rad);
        L_k_est = min(L_k_est, 100);  % 上限100m（物理合理性约束）
    else
        L_k_est = 0;
    end

    fprintf('  Angle width (-6dB): %.2f°, L_k est: %.3f m\n', angle_width_deg, L_k_est);

    %% --- 3e. 散射中心分类 ---
    % 判别树 (来自md文档)
    isSliding = (delta_R_norm > slidingThreshold) && (R2_sliding > 0.5);

    % 距离展宽检查
    hr_profile = residual_HRRP_mag(peak_ang_idx, :);  % 峰值角度处的距离像
    hr_dB = 20 * log10(hr_profile / max(hr_profile) + epsilon_dB);
    above_6dB_range = hr_dB > -6;
    range_width_bins = sum(above_6dB_range);
    range_width_m = range_width_bins * delta_r_fft;
    range_width_norm = range_width_m / delta_r;

    % 角度展宽检查
    if angle_width_deg > 0
        % -6dB角度宽度 vs 分辨单元对应的sinc宽度
        % sinc宽度 = c/(f_c * D) 其中D是目标尺寸
        angle_width_norm = angle_width_deg / 2;  % 简化判断
    else
        angle_width_norm = 0;
    end

    if isSliding
        type_str = 'sliding';
        type_cn  = '滑动型';
    elseif range_width_norm > distributedWidthFactor
        type_str = 'distributed';
        type_cn  = '分布型';
    elseif angle_width_deg > 2 * delt  % 角度展宽 > 2×角度步长
        type_str = 'distributed';
        type_cn  = '分布型';
    else
        type_str = 'local';
        type_cn  = '局部型';
    end

    % α值物理含义对照
    alpha_mechanism = getAlphaMechanism(alpha_k);

    fprintf('  Classification: %s (%s)\n', type_cn, type_str);
    fprintf('  α=%.2f → %s\n', alpha_k, alpha_mechanism);
    fprintf('  Range width: %.1f×Δr, Angle width: %.1f°\n', range_width_norm, angle_width_deg);

    %% --- 3f. 构建散射中心回波并减去 ---
    phi_bar = phi_k;  % 中心角度
    L_k = L_k_est;

    % 合成回波 S_k(f, φ) for all angles
    S_k = synthesizeScattererEcho(freq_array, angle_vec, f_c, c, ...
        A_k, alpha_k, R_k, L_k, phi_bar, dR_dphi, type_str, cleanGain, peak_ang_idx);

    % 从残差中减去
    residual = residual - S_k;

    %% --- 3g. 保存散射中心参数 ---
    K = K + 1;
    centers(K).A        = A_k;        %#ok<AGROW>
    centers(K).alpha    = alpha_k;    %#ok<AGROW>
    centers(K).R        = R_k;        %#ok<AGROW>
    centers(K).L        = L_k;        %#ok<AGROW>
    centers(K).phi_bar  = phi_bar;    %#ok<AGROW>
    centers(K).dR_dphi  = dR_dphi;    %#ok<AGROW>
    centers(K).type     = type_str;   %#ok<AGROW>
    centers(K).type_cn  = type_cn;    %#ok<AGROW>
    centers(K).alpha_mechanism = alpha_mechanism; %#ok<AGROW>
    centers(K).angle_idx = peak_ang_idx; %#ok<AGROW>
    centers(K).range_idx = peak_r_idx;   %#ok<AGROW>
    centers(K).R2_alpha = R2;         %#ok<AGROW>
    centers(K).R2_sliding = R2_sliding; %#ok<AGROW>
    centers(K).angle_width_deg = angle_width_deg; %#ok<AGROW>
    centers(K).range_width_m = range_width_m;     %#ok<AGROW>

    clean_history(iter).K = K;         %#ok<AGROW>
    clean_history(iter).residual_energy = residual_energy; %#ok<AGROW>
    clean_history(iter).energy_ratio = energy_ratio; %#ok<AGROW>
end

fprintf('\n  CLEAN finished: %d scattering centers extracted\n', K);

if K == 0
    error('No scattering centers extracted. Try lowering CleanThreshold.');
end

%% ========================================================================
%% 步骤5: 参数优化
%% ========================================================================
fprintf('\n========================================\n');
fprintf('  Step 5: Parameter Optimization\n');
fprintf('  (参数优化)\n');
fprintf('========================================\n');

if doOptimize && K > 0
    fprintf('  Optimizing %d scattering centers jointly...\n', K);

    % 构建初始参数向量
    x0 = centersToVector(centers);

    % 参数边界和约束
    [lb, ub] = getParamBounds(centers, f_c, delta_r, angle_vec);

    % Levenberg-Marquardt 非线性最小二乘优化
    opts = optimoptions('lsqnonlin', ...
        'Display', 'iter-detailed', ...
        'MaxIterations', 200, ...
        'MaxFunctionEvaluations', 5000, ...
        'OptimalityTolerance', 1e-8, ...
        'StepTolerance', 1e-10, ...
        'FunctionTolerance', 1e-8, ...
        'Algorithm', 'levenberg-marquardt');

    cost_func = @(x) computeResidual(x, S_complex, freq_array, angle_vec, ...
        f_c, c, K, centers);

    try
        [x_opt, resnorm, residual_vec, exitflag] = ...
            lsqnonlin(cost_func, x0, lb, ub, opts);

        centers = vectorToCenters(x_opt, centers);

        fprintf('  Optimization complete: exitflag=%d, residual norm=%.4e\n', ...
            exitflag, resnorm);

        % 更新合成回波
        S_model_opt = synthesizeFullEcho(freq_array, angle_vec, f_c, c, centers);
        residual_opt = S_complex - S_model_opt;
        final_residual_ratio = sum(abs(residual_opt(:)).^2) / initial_energy;
        fprintf('  Final residual energy: %.2f%% of initial\n', 100*final_residual_ratio);

    catch ME
        warning('Optimization failed: %s. Using initial parameters.', ME.message);
        S_model_opt = synthesizeFullEcho(freq_array, angle_vec, f_c, c, centers);
        residual_opt = S_complex - S_model_opt;
    end
else
    fprintf('  Skipping optimization.\n');
    S_model_opt = synthesizeFullEcho(freq_array, angle_vec, f_c, c, centers);
    residual_opt = S_complex - S_model_opt;
end

%% ========================================================================
%% 步骤6: 模型验证
%% ========================================================================
fprintf('\n========================================\n');
fprintf('  Step 6: Model Validation\n');
fprintf('  (模型验证)\n');
fprintf('========================================\n');

% 计算HRRP用于验证
HRRP_po_dB = HRRP_dB;  % 原始PO-HRRP
HRRP_model_complex = zeros(N_angles, N_fft);
for i_ang = 1:N_angles
    Es_f = S_model_opt(i_ang, :);
    Es_windowed = W_freq .* Es_f;
    Es_padded = [Es_windowed, zeros(1, N_fft - N_f)];
    h = ifft(Es_padded, N_fft);
    h = ifftshift(h);
    HRRP_model_complex(i_ang, :) = h;
end
HRRP_model_dB = 20 * log10(abs(HRRP_model_complex) + epsilon_dB);

% 全局归一化（PO和模型使用相同的参考）
global_peak = max(HRRP_po_dB(:));
HRRP_po_disp = HRRP_po_dB - global_peak;
HRRP_model_disp = HRRP_model_dB - global_peak;

%% --- 6a. 回波幅度相关系数 ---
S_po_vec = S_complex(:);
S_model_vec = S_model_opt(:);
rho = abs(S_po_vec' * S_model_vec) / ...
    sqrt((S_po_vec' * S_po_vec) * (S_model_vec' * S_model_vec));

fprintf('  Correlation coefficient ρ = %.4f\n', rho);
if rho > 0.9
    fprintf('  → Excellent agreement (ρ > 0.9)\n');
elseif rho > 0.8
    fprintf('  → Good agreement (ρ > 0.8)\n');
elseif rho > 0.7
    fprintf('  → Fair agreement (ρ > 0.7)\n');
else
    fprintf('  → Poor agreement (ρ < 0.7), consider more centers or wider angle range\n');
end

%% --- 6b. 逐角度RCS对比 ---
PO_RCS_dB = 20 * log10(sqrt(mean(abs(S_complex).^2, 2)) + epsilon_dB);
Model_RCS_dB = 20 * log10(sqrt(mean(abs(S_model_opt).^2, 2)) + epsilon_dB);

%% --- 6c. 散射中心参数表 ---
fprintf('\n  ========================================\n');
fprintf('  Scattering Center Parameter Table\n');
fprintf('  (散射中心参数表)\n');
fprintf('  ========================================\n');
fprintf('  %-4s %-10s %-10s %-10s %-10s %-12s %-12s\n', ...
    'ID', 'Type', 'A(dB)', 'α', 'R(m)', 'L(m)', 'dR/dφ(m/rad)', 'φ̄(°)');
fprintf('  %s\n', repmat('-', 1, 85));
for k = 1:K
    fprintf('  %-4d %-10s %-10.1f %-10.3f %-10.4f %-10.4f %-12.4f %-12.1f\n', ...
        k, centers(k).type, 20*log10(centers(k).A), centers(k).alpha, ...
        centers(k).R, centers(k).L, centers(k).dR_dphi, centers(k).phi_bar);
end

%% ========================================================================
%% 绘图
%% ========================================================================
fprintf('\n========================================\n');
fprintf('  Generating Validation Plots\n');
fprintf('  (生成验证图)\n');
fprintf('========================================\n');

dyn_range = 60;  % 显示动态范围 (dB)

% 动态范围裁剪
po_disp = HRRP_po_disp;
model_disp = HRRP_model_disp;
po_disp(po_disp < -dyn_range) = -dyn_range;
model_disp(model_disp < -dyn_range) = -dyn_range;

%% --- 图1: HRRP历程对比（四联图）---
figure(1);
clf;
set(gcf, 'Name', 'HRRP Map Comparison (历程图对比)', ...
    'NumberTitle', 'off', 'Position', [50, 80, 1400, 900]);

% (a) PO原始HRRP历程图
subplot(2, 3, 1);
imagesc(angle_vec, range_axis, po_disp');
colormap('jet');
axis xy;
xlabel(sprintf('Angle %s (%s)', angleLabel, angleUnit), 'FontSize', 11);
ylabel('Range (m)', 'FontSize', 11);
title(sprintf('PO HRRP Map\n(PO一维距离像历程图)'), 'FontSize', 12, 'FontWeight', 'bold');
clim([-dyn_range, 0]);
colorbar;
set(gca, 'FontSize', 10);

% (b) 散射中心模型HRRP历程图
subplot(2, 3, 2);
imagesc(angle_vec, range_axis, model_disp');
colormap('jet');
axis xy;
xlabel(sprintf('Angle %s (%s)', angleLabel, angleUnit), 'FontSize', 11);
ylabel('Range (m)', 'FontSize', 11);
title(sprintf('Model HRRP Map (%d centers)\n(模型一维距离像历程图)', K), ...
    'FontSize', 12, 'FontWeight', 'bold');
clim([-dyn_range, 0]);
colorbar;
set(gca, 'FontSize', 10);

% (c) 残差HRRP历程图
residual_HRRP_dB = 20 * log10(abs(HRRP_complex - HRRP_model_complex) + epsilon_dB);
residual_disp = residual_HRRP_dB - global_peak;
residual_disp(residual_disp < -dyn_range) = -dyn_range;

subplot(2, 3, 3);
imagesc(angle_vec, range_axis, residual_disp');
colormap('jet');
axis xy;
xlabel(sprintf('Angle %s (%s)', angleLabel, angleUnit), 'FontSize', 11);
ylabel('Range (m)', 'FontSize', 11);
title(sprintf('Residual HRRP Map (ρ=%.3f)\n(残差历程图)', rho), ...
    'FontSize', 12, 'FontWeight', 'bold');
clim([-dyn_range, 0]);
colorbar;
set(gca, 'FontSize', 10);

% (d) 逐角度RCS对比
subplot(2, 3, 4);
plot(angle_vec, PO_RCS_dB, 'b-', 'LineWidth', 1.5);
hold on;
plot(angle_vec, Model_RCS_dB, 'r--', 'LineWidth', 1.5);
hold off;
xlabel(sprintf('Angle %s (%s)', angleLabel, angleUnit), 'FontSize', 11);
ylabel('RCS (dBsm)', 'FontSize', 11);
title('Angle-by-Angle RCS Comparison (逐角度RCS对比)', 'FontSize', 12, 'FontWeight', 'bold');
legend('PO', sprintf('Model (%d centers)', K), 'Location', 'best');
grid on;
set(gca, 'FontSize', 10);

% (e) 单角度距离像对比（取中间角度）
mid_angle_idx = round(N_angles / 2);
subplot(2, 3, 5);
plot(range_axis, HRRP_po_disp(mid_angle_idx, :), 'b-', 'LineWidth', 1.2);
hold on;
plot(range_axis, HRRP_model_disp(mid_angle_idx, :), 'r--', 'LineWidth', 1.2);
% 标注散射中心位置
for k = 1:K
    if abs(angle_vec(mid_angle_idx) - centers(k).phi_bar) < 30
        xline(centers(k).R, 'g--', sprintf('#%d', k), 'LineWidth', 1, 'Alpha', 0.7);
    end
end
hold off;
xlabel('Range (m)', 'FontSize', 11);
ylabel('Relative Intensity (dB)', 'FontSize', 11);
title(sprintf('1D Range Profile at %s=%.1f°\n(单角度一维距离像)', ...
    angleLabel, angle_vec(mid_angle_idx)), 'FontSize', 12, 'FontWeight', 'bold');
legend('PO', 'Model', 'Location', 'best');
grid on;
set(gca, 'FontSize', 10);

% (f) CLEAN收敛曲线
subplot(2, 3, 6);
if ~isempty(clean_history)
    energy_ratios = [clean_history.energy_ratio] * 100;
    semilogy(1:length(energy_ratios), energy_ratios, 'bo-', ...
        'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    xlabel('CLEAN Iteration', 'FontSize', 11);
    ylabel('Residual Energy (%)', 'FontSize', 11);
    title('CLEAN Convergence (CLEAN收敛曲线)', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 10);
    yline(cleanThreshold*100, 'r--', 'Threshold', 'LineWidth', 1);
end

sgtitle(sprintf(['Scattering Center Modeling: %s\n', ...
    '%d Centers, ρ=%.3f, BW=%.1f GHz, Δr=%.2f cm'], ...
    modelName, K, rho, B/1e9, delta_r*100), ...
    'FontSize', 14, 'FontWeight', 'bold');

%% --- 图2: 散射中心参数可视化 ---
if showDetailed
    figure(2);
    clf;
    set(gcf, 'Name', 'Scattering Center Parameters (散射中心参数)', ...
        'NumberTitle', 'off', 'Position', [100, 100, 1200, 600]);

    % (a) 散射中心在距离-角度平面上的分布
    subplot(2, 3, 1);
    type_colors = containers.Map();
    type_colors('local')       = [0.0, 0.45, 0.74];  % 蓝
    type_colors('distributed') = [0.85, 0.33, 0.10];  % 橙
    type_colors('sliding')     = [0.93, 0.69, 0.13];  % 黄

    hold on;
    for k = 1:K
        c = type_colors(centers(k).type);
        h = plot(centers(k).phi_bar, centers(k).R, 'o', ...
            'MarkerSize', 8 + 12 * centers(k).A / max([centers.A]), ...
            'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    end
    hold off;
    xlabel(sprintf('Angle %s (%s)', angleLabel, angleUnit), 'FontSize', 11);
    ylabel('Range R (m)', 'FontSize', 11);
    title('Scatterers in Angle-Range Plane', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 10);
    % 图例
    legend_str = {};
    for t = {'local', 'distributed', 'sliding'}
        if any(strcmp({centers.type}, t{1}))
            legend_str{end+1} = t{1}; %#ok<AGROW>
        end
    end
    legend(legend_str, 'Location', 'best');

    % (b) α值分布（按散射中心编号）
    subplot(2, 3, 2);
    alpha_vals = [centers.alpha];
    bar(1:K, alpha_vals, 'FaceColor', [0.3, 0.6, 0.9]);
    % 标注物理含义
    hold on;
    yline(1, 'r--', 'α=1 (尖端)', 'LineWidth', 1, 'Alpha', 0.5);
    yline(0.5, 'g--', 'α=0.5 (边缘)', 'LineWidth', 1, 'Alpha', 0.5);
    yline(0, 'm--', 'α=0 (镜面)', 'LineWidth', 1, 'Alpha', 0.5);
    yline(-0.5, 'c--', 'α=-0.5 (长边)', 'LineWidth', 1, 'Alpha', 0.5);
    yline(-1, 'k--', 'α=-1 (角顶)', 'LineWidth', 1, 'Alpha', 0.5);
    hold off;
    xlabel('Scatterer ID', 'FontSize', 11);
    ylabel('Frequency Dependence α', 'FontSize', 11);
    title('α Values (频率依赖指数)', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 10);

    % (c) 散射中心类型饼图
    subplot(2, 3, 3);
    type_counts = [sum(strcmp({centers.type}, 'local')), ...
                   sum(strcmp({centers.type}, 'distributed')), ...
                   sum(strcmp({centers.type}, 'sliding'))];
    if all(type_counts == 0)
        text(0.5, 0.5, 'No centers', 'HorizontalAlignment', 'center');
    else
        pie(type_counts, {'Local (局部型)', 'Distributed (分布型)', 'Sliding (滑动型)'});
        colormap(gca, [type_colors('local'); type_colors('distributed'); type_colors('sliding')]);
    end
    title('Scatterer Type Distribution', 'FontSize', 12, 'FontWeight', 'bold');

    % (d) 幅度分布
    subplot(2, 3, 4);
    A_dB_vals = 20 * log10([centers.A]);
    bar(1:K, A_dB_vals, 'FaceColor', [0.9, 0.5, 0.2]);
    xlabel('Scatterer ID', 'FontSize', 11);
    ylabel('Amplitude (dB)', 'FontSize', 11);
    title('Scatterer Amplitudes (幅度分布)', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    set(gca, 'FontSize', 10);

    % (e) 散射中心物理含义速查
    subplot(2, 3, 5);
    axis off;
    text(0.05, 0.95, 'α 物理含义速查:', 'FontSize', 12, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'top');
    mechanism_table = {
        'α = 1.0   尖端绕射 / 双曲率曲面 (头锥尖端)';
        'α = 0.5   边缘绕射 / 单曲率曲面 (锥柱连接处)';
        'α = 0.0   平板镜面反射 / 角反射器 (台体端面)';
        'α = -0.5  长直边缘掠入射 (柱体边缘)';
        'α = -1.0  角顶绕射 (非常尖锐的角)';
    };
    for i = 1:length(mechanism_table)
        text(0.05, 0.85 - i*0.15, mechanism_table{i}, 'FontSize', 9, ...
            'VerticalAlignment', 'top');
    end

    % 模型统计信息
    text(0.05, 0.05, sprintf(['Statistics:\n', ...
        '  Total centers: %d\n', ...
        '  Local: %d, Distributed: %d, Sliding: %d\n', ...
        '  ρ = %.4f\n', ...
        '  α range: [%.2f, %.2f]'], ...
        K, type_counts(1), type_counts(2), type_counts(3), ...
        rho, min(alpha_vals), max(alpha_vals)), ...
        'FontSize', 9, 'VerticalAlignment', 'bottom');

    % (f) 残差能量 vs 散射中心数量
    subplot(2, 3, 6);
    if ~isempty(clean_history)
        energies = [clean_history.residual_energy] / initial_energy * 100;
        plot(1:length(energies), energies, 'bs-', ...
            'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
        xlabel('Number of Extracted Centers', 'FontSize', 11);
        ylabel('Residual Energy Ratio (%)', 'FontSize', 11);
        title('Residual vs. Number of Centers', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        set(gca, 'FontSize', 10);
    end

    sgtitle(sprintf('Scattering Center Analysis: %s (%d centers, ρ=%.3f)', ...
        modelName, K, rho), 'FontSize', 14, 'FontWeight', 'bold');
end

%% --- 图3: 频域回波对比 ---
if showDetailed
    figure(3);
    clf;
    set(gcf, 'Name', 'Frequency-Domain Echo Comparison (频域回波对比)', ...
        'NumberTitle', 'off', 'Position', [150, 150, 1000, 600]);

    % 选取3个代表性角度
    angles_to_show = [1, round(N_angles/2), N_angles];
    for i_show = 1:3
        i_ang = angles_to_show(i_show);

        subplot(3, 2, 2*i_show-1);
        S_po_db = 20 * log10(abs(S_complex(i_ang, :)) + epsilon_dB);
        S_model_db = 20 * log10(abs(S_model_opt(i_ang, :)) + epsilon_dB);
        plot(freq_array/1e9, S_po_db, 'b-', 'LineWidth', 1.2);
        hold on;
        plot(freq_array/1e9, S_model_db, 'r--', 'LineWidth', 1.2);
        hold off;
        xlabel('Frequency (GHz)', 'FontSize', 10);
        ylabel('|E_s| (dB)', 'FontSize', 10);
        title(sprintf('Spectrum at %s=%.1f°', angleLabel, angle_vec(i_ang)), ...
            'FontSize', 11, 'FontWeight', 'bold');
        legend('PO', 'Model', 'Location', 'best');
        grid on;
        set(gca, 'FontSize', 9);

        subplot(3, 2, 2*i_show);
        phase_po = unwrap(angle(S_complex(i_ang, :)));
        phase_model = unwrap(angle(S_model_opt(i_ang, :)));
        plot(freq_array/1e9, phase_po, 'b-', 'LineWidth', 1.2);
        hold on;
        plot(freq_array/1e9, phase_model, 'r--', 'LineWidth', 1.2);
        hold off;
        xlabel('Frequency (GHz)', 'FontSize', 10);
        ylabel('Phase (rad)', 'FontSize', 10);
        title(sprintf('Phase at %s=%.1f°', angleLabel, angle_vec(i_ang)), ...
            'FontSize', 11, 'FontWeight', 'bold');
        legend('PO', 'Model', 'Location', 'best');
        grid on;
        set(gca, 'FontSize', 9);
    end

    sgtitle(sprintf('Frequency-Domain Validation: %s (ρ=%.3f)', modelName, rho), ...
        'FontSize', 14, 'FontWeight', 'bold');
end

%% ========================================================================
%% 保存结果
%% ========================================================================
fprintf('\nSaving results...\n');

% 创建输出目录
[resultDir, nowStr] = createResultDir('main_scattering_center_modeling');

% --- 保存图像 ---
figFile1 = fullfile(resultDir, ['hrrp_validation_' nowStr '.png']);
saveas(1, figFile1);
fprintf('  Figure 1 (HRRP validation): %s\n', figFile1);

if showDetailed
    figFile2 = fullfile(resultDir, ['scatterer_params_' nowStr '.png']);
    saveas(2, figFile2);
    fprintf('  Figure 2 (Parameters):      %s\n', figFile2);

    figFile3 = fullfile(resultDir, ['freq_domain_validation_' nowStr '.png']);
    saveas(3, figFile3);
    fprintf('  Figure 3 (Freq domain):     %s\n', figFile3);
end

% --- 保存散射中心参数表为文本文件 ---
paramFile = fullfile(resultDir, ['scattering_centers_' nowStr '.txt']);
fid = fopen(paramFile, 'w');
fprintf(fid, '# Scattering Center Model Parameters\n');
fprintf(fid, '# Generated: %s\n', nowStr);
fprintf(fid, '# Model: %s\n', modelName);
fprintf(fid, '# Frequency: %.2f - %.2f GHz, BW: %.1f GHz\n', ...
    freq_array(1)/1e9, freq_array(end)/1e9, B/1e9);
fprintf(fid, '# Range resolution: %.4f m\n', delta_r);
fprintf(fid, '# Correlation coefficient ρ = %.4f\n', rho);
fprintf(fid, '# Number of scattering centers: %d\n\n', K);
fprintf(fid, '# Columns: ID  Type  A(dB)  Alpha  R(m)  L(m)  dR_dphi(m/rad)  phi_bar(deg)  Alpha_Mechanism\n');
fprintf(fid, '# %s\n', repmat('-', 1, 95));
for k = 1:K
    fprintf(fid, '%-4d  %-12s  %8.2f  %7.3f  %9.4f  %9.4f  %14.4f  %12.1f  %s\n', ...
        k, centers(k).type, 20*log10(centers(k).A), centers(k).alpha, ...
        centers(k).R, centers(k).L, centers(k).dR_dphi, centers(k).phi_bar, ...
        centers(k).alpha_mechanism);
end
fclose(fid);
fprintf('  Parameter file:             %s\n', paramFile);

% --- 保存 .mat 数据文件 ---
matFile = fullfile(resultDir, ['scattering_centers_' nowStr '.mat']);
save(matFile, 'centers', 'K', 'rho', 'S_model_opt', 'residual_opt', ...
    'freq_array', 'angle_vec', 'f_c', 'c', 'B', 'delta_r', 'angleLabel', ...
    'modelName', 'clean_history', 'initial_energy', '-v7.3');
fprintf('  Data file (.mat):           %s\n', matFile);

%% ========================================================================
%% 完成摘要
%% ========================================================================
fprintf('\n========================================\n');
fprintf('  Scattering Center Modeling Complete\n');
fprintf('  (散射中心建模完成)\n');
fprintf('========================================\n');
fprintf('  Model:             %s\n', modelName);
fprintf('  Bandwidth:         %.2f GHz\n', B/1e9);
fprintf('  Range resolution:  %.4f m\n', delta_r);
fprintf('  Scattering centers: %d\n', K);
fprintf('    - Local:         %d\n', sum(strcmp({centers.type}, 'local')));
fprintf('    - Distributed:   %d\n', sum(strcmp({centers.type}, 'distributed')));
fprintf('    - Sliding:       %d\n', sum(strcmp({centers.type}, 'sliding')));
fprintf('  Correlation ρ:    %.4f\n', rho);
fprintf('  Output directory:  %s\n', resultDir);
fprintf('========================================\n');

end  % main_scattering_center_modeling

%% ========================================================================
%% 子函数: 合成单个散射中心的回波
%% ========================================================================
function S_k = synthesizeScattererEcho(freq_array, angle_vec, f_c, c, ...
    A, alpha, R, L, phi_bar, dR_dphi, type, gain, peak_angle_idx)
% 合成单个散射中心的频域回波 S_k(f, φ)
%
% Input:
%   freq_array    - 频率数组 (1 × N_f) Hz
%   angle_vec     - 角度数组 (N_angles × 1)
%   f_c           - 中心频率 (Hz)
%   c             - 光速 (m/s)
%   A             - 幅度
%   alpha         - 频率依赖指数
%   R             - 距离位置 (m)
%   L             - 有效长度 (m, 分布型)
%   phi_bar       - 中心角度 (deg)
%   dR_dphi       - 距离滑动率 (m/rad, 滑动型)
%   type          - 类型: 'local', 'distributed', 'sliding'
%   gain          - CLEAN增益因子
%   peak_angle_idx - 峰值所在的角度索引
%
% Output:
%   S_k           - 合成回波 (N_angles × N_f)

N_f = length(freq_array);
N_angles = length(angle_vec);
phi_rad = deg2rad(angle_vec);
phi_bar_rad = deg2rad(phi_bar);

% 频率因子 (j*f/f_c)^α
freq_factor = (1j * freq_array / f_c).^alpha;

% 初始化
S_k = zeros(N_angles, N_f);

switch type
    case 'local'
        % 局部型: 窄角度响应 (高斯窗口)
        % 角度宽度自适应：根据实际峰值宽度
        sigma_angle = deg2rad(max(2, 180/N_angles * 2));  % 至少覆盖2个角度采样点
        angle_window = exp(-(phi_rad - phi_bar_rad).^2 / (2 * sigma_angle^2));

        for i_ang = 1:N_angles
            if angle_window(i_ang) < 0.01
                continue;  % 跳过贡献极小的角度
            end
            phase_term = exp(-1j * (4*pi/c) * R * freq_array);
            S_k(i_ang, :) = gain * A * angle_window(i_ang) * freq_factor .* phase_term;
        end

    case 'distributed'
        % 分布型: sinc 角度响应
        for i_ang = 1:N_angles
            dphi = phi_rad(i_ang) - phi_bar_rad;
            if abs(dphi) < 1e-10
                sinc_val = 1;
            else
                arg = (2*pi/c) * freq_array * L * sin(dphi);
                % 向量化的sinc
                sinc_val = ones(size(arg));
                nonzero = abs(arg) > 1e-15;
                sinc_val(nonzero) = sin(arg(nonzero)) ./ arg(nonzero);
            end

            phase_term = exp(-1j * (4*pi/c) * R * freq_array);
            S_k(i_ang, :) = gain * A * freq_factor .* phase_term .* sinc_val;
        end

    case 'sliding'
        % 滑动型: 距离位置随角度变化 R(φ) = R + dR/dφ * (φ - φ̄)
        % 窄角度窗口
        sigma_angle = deg2rad(max(2, 180/N_angles * 2));
        angle_window = exp(-(phi_rad - phi_bar_rad).^2 / (2 * sigma_angle^2));

        for i_ang = 1:N_angles
            if angle_window(i_ang) < 0.01
                continue;
            end
            R_phi = R + dR_dphi * (phi_rad(i_ang) - phi_bar_rad);
            phase_term = exp(-1j * (4*pi/c) * R_phi * freq_array);
            S_k(i_ang, :) = gain * A * angle_window(i_ang) * freq_factor .* phase_term;
        end

    otherwise
        error('Unknown scatterer type: %s', type);
end
end

%% ========================================================================
%% 子函数: 合成所有散射中心的总回波
%% ========================================================================
function S_total = synthesizeFullEcho(freq_array, angle_vec, f_c, c, centers)
% 合成所有散射中心的总回波
N_f = length(freq_array);
N_angles = length(angle_vec);
K = length(centers);

S_total = zeros(N_angles, N_f);

for k = 1:K
    ct = centers(k);
    S_k = synthesizeScattererEcho(freq_array, angle_vec, f_c, c, ...
        ct.A, ct.alpha, ct.R, ct.L, ct.phi_bar, ct.dR_dphi, ct.type, 1.0, ct.angle_idx);
    S_total = S_total + S_k;
end
end

%% ========================================================================
%% 子函数: 参数向量 ↔ 结构体转换
%% ========================================================================
function x = centersToVector(centers)
% 将散射中心结构体数组转换为优化参数向量
% 每个散射中心: [A, alpha, R, L, dR_dphi]
K = length(centers);
x = zeros(5 * K, 1);
for k = 1:K
    idx = (k-1)*5 + 1;
    x(idx)   = centers(k).A;
    x(idx+1) = centers(k).alpha;
    x(idx+2) = centers(k).R;
    x(idx+3) = centers(k).L;
    x(idx+4) = centers(k).dR_dphi;
end
end

function centers = vectorToCenters(x, centers_template)
% 将优化参数向量转换回散射中心结构体数组
K = length(centers_template);
centers = centers_template;
for k = 1:K
    idx = (k-1)*5 + 1;
    centers(k).A       = x(idx);
    centers(k).alpha   = x(idx+1);
    centers(k).R       = x(idx+2);
    centers(k).L       = x(idx+3);
    centers(k).dR_dphi = x(idx+4);
end
end

%% ========================================================================
%% 子函数: 参数边界
%% ========================================================================
function [lb, ub] = getParamBounds(centers, f_c, delta_r, angle_vec)
% 获取参数优化的边界
K = length(centers);
lb = zeros(5 * K, 1);
ub = zeros(5 * K, 1);

for k = 1:K
    idx = (k-1)*5 + 1;

    % A: 幅度 (正值)
    A0 = centers(k).A;
    lb(idx)   = A0 * 0.1;
    ub(idx)   = A0 * 10;

    % alpha: 频率依赖指数 [-1.5, 1.5]
    lb(idx+1) = -1.5;
    ub(idx+1) = 1.5;

    % R: 距离位置 (在当前位置附近 ±5个分辨单元)
    R0 = centers(k).R;
    lb(idx+2) = R0 - 5 * delta_r;
    ub(idx+2) = R0 + 5 * delta_r;

    % L: 有效长度 [0, 100]
    lb(idx+3) = 0;
    ub(idx+3) = 100;

    % dR_dphi: 滑动率 [-10, 10]
    lb(idx+4) = -10;
    ub(idx+4) = 10;
end
end

%% ========================================================================
%% 子函数: 优化残差计算
%% ========================================================================
function r = computeResidual(x, S_measured, freq_array, angle_vec, ...
    f_c, c, K, centers_template)
% 计算模型合成回波与测量回波之间的残差向量
centers_current = vectorToCenters(x, centers_template);

S_model = synthesizeFullEcho(freq_array, angle_vec, f_c, c, centers_current);

% 复数残差（实部和虚部分开，因为lsqnonlin需要实数）
diff = S_model - S_measured;
r = [real(diff(:)); imag(diff(:))];
end

%% ========================================================================
%% 子函数: α值物理含义
%% ========================================================================
function mechanism = getAlphaMechanism(alpha)
% 根据α值判断最可能的散射机理
if alpha > 0.75
    mechanism = '尖端绕射/双曲率曲面';
elseif alpha > 0.25
    mechanism = '边缘绕射/单曲率曲面';
elseif alpha > -0.25
    mechanism = '镜面反射/角反射器';
elseif alpha > -0.75
    mechanism = '长直边缘(掠入射)';
else
    mechanism = '角顶绕射';
end
end

%% ========================================================================
%% 子函数: 三元运算符辅助函数
%% ========================================================================
function result = ternary(condition, true_val, false_val)
if condition
    result = true_val;
else
    result = false_val;
end
end
