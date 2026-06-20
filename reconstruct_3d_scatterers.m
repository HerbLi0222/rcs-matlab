function reconstruct_3d_scatterers(varargin)
% RECONSTRUCT_3D_SCATTERERS  三维散射中心重建与可视化
%
%   从 wideband_scattering 数据中，在多个中心观测角度下进行 SAR 成像和
%   CLEAN 散射中心提取，将 (range, cross-range) 二维坐标逆投影回三维空间，
%   并与原始目标 STL 模型叠加显示。
%
%   方法原理:
%     散射仿真使用 exp(+j·2k·R·P) 相位约定（标准雷达用 exp(-j·2k·R·P)），
%     导致 MATLAB IFFT 后 (r, x) 坐标发生符号翻转：
%       r_sar = -k̂·P,   x_sar = -k̂_perp·P
%
%     在中心角 θ₀ (phi=0) 的 SAR 图像中：
%       k̂       = [sin(θ₀), 0, cos(θ₀)]   — 雷达视线方向 (Range)
%       k̂_perp  = [cos(θ₀), 0, -sin(θ₀)]  — 方位向 (Cross-Range)
%
%     修正后的逆变换（从2D到3D）：
%       P_3D = - r_sar·k̂ - x_sar·k̂_perp
%     （相位参考为原点，重建结果直接位于 STL 模型坐标系）
%
%     注意：单平面扫描无法分辨垂直于扫描平面（y轴）的分量，
%     因此重建的散射点始终位于 y=0 平面（xz扫描平面）内。
%
%   Usage:
%     >> reconstruct_3d_scatterers
%     >> reconstruct_3d_scatterers('DataFile', 'wideband_scattering_xxx.mat')
%     >> reconstruct_3d_scatterers('CenterAngles', [30, 60, 90, 120, 150])
%     >> reconstruct_3d_scatterers('CleanMaxIter', 30, 'CleanThreshold', 0.03)

%% ========================================================================
%% 0. 参数解析
%% ========================================================================
addpath('lib');

paramNames = {'DataFile','CenterAngles','ApertureDeg','CleanMaxIter',...
    'CleanThreshold','CleanGain','ZeroPadFactor','WindowType'};

p = inputParser;
p.addOptional('DataFile', '', @(x) ischar(x) || isstring(x) || isempty(x));
p.addParameter('CenterAngles', 0:15:360, @(x) isnumeric(x));
p.addParameter('ApertureDeg', 10, @(x) isnumeric(x) && isscalar(x));
p.addParameter('CleanMaxIter', 30, @(x) isnumeric(x) && isscalar(x));
p.addParameter('CleanThreshold', 0.05, @(x) isnumeric(x) && isscalar(x));
p.addParameter('CleanGain', 0.6, @(x) isnumeric(x) && isscalar(x));
p.addParameter('ZeroPadFactor', 2, @(x) isnumeric(x) && isscalar(x));
p.addParameter('WindowType', 'chebwin', @(x) ischar(x) || isstring(x));

if ~isempty(varargin) && (ischar(varargin{1}) || isstring(varargin{1})) ...
        && any(strcmp(char(varargin{1}), paramNames))
    varargin = [{''}, varargin];
end
p.parse(varargin{:});

resultFile       = p.Results.DataFile;
center_angles    = p.Results.CenterAngles;
aperture_deg     = p.Results.ApertureDeg;
cleanMaxIter     = p.Results.CleanMaxIter;
cleanThreshold   = p.Results.CleanThreshold;
cleanGain        = p.Results.CleanGain;
zeroPadFactor    = p.Results.ZeroPadFactor;
windowType       = p.Results.WindowType;

%% ========================================================================
%% 1. 加载宽带散射数据和目标几何
%% ========================================================================
fprintf('========================================\n');
fprintf('  3D Scattering Center Reconstruction\n');
fprintf('  (三维散射中心重建)\n');
fprintf('========================================\n\n');

% --- 1a. 加载 wideband_scattering 数据 ---
if isempty(resultFile)
    resultFile = findLatestResultFile('wideband_scattering_*.mat');
    if isempty(resultFile)
        error('No wideband_scattering_*.mat found in results/.\n');
    end
end

fprintf('Loading wideband scattering data...\n');
fprintf('  File: %s\n', resultFile);
data = load(resultFile);

S_complex      = data.S_complex;
freq_array     = data.freq_array(:)';
theta_array    = data.theta_array;
phi_array      = data.phi_array;
c              = data.c;
B              = data.B;
delta_r        = data.delta_r;
ip             = data.ip;
it             = data.it;
N_f            = data.N_f;

% 中心频率和波长
f_center = (freq_array(1) + freq_array(end)) / 2;
lambda_c = c / f_center;

% 模型名和包围盒中心
if isfield(data, 'inputModel')
    [~, modelName, ~] = fileparts(data.inputModel);
    modelFile = fullfile('stl_models', data.inputModel);
else
    modelName = 'Unknown';
    modelFile = '';
end

if isfield(data, 'bbox_center')
    P_ref = data.bbox_center(:);  % 3x1
    fprintf('  Phase reference (bbox center): [%.3f, %.3f, %.3f] m\n', P_ref);
else
    P_ref = [0; 0; 0];
    fprintf('  Phase reference: origin (no bbox data)\n');
end

% 角度向量
if isfield(data, 'delt'), delt = data.delt; else delt = 1; end
if isfield(data, 'tstart'), tstart = data.tstart; else tstart = 0; end
if ip == 1 && it > 1
    theta_vec = theta_array(1, :)';
elseif it == 1 && ip > 1
    theta_vec = phi_array(:, 1);
else
    theta_vec = (1:size(S_complex, 1))';
end

fprintf('  Model: %s, Bandwidth: %.1f GHz, λc: %.4f m\n', ...
    modelName, B/1e9, lambda_c);
fprintf('  Freq points: %d, Angle range: [%.0f, %.0f]°\n', ...
    N_f, theta_vec(1), theta_vec(end));

% --- 1b. 加载 STL 模型几何 ---
fprintf('\nLoading target geometry...\n');
[vertices, faces] = readSTLGeometry(modelFile, data);
fprintf('  Vertices: %d, Faces: %d\n', size(vertices, 1), size(faces, 1));

%% ========================================================================
%% 2. 多角度 SAR 成像与散射中心提取
%% ========================================================================
N_angles = length(center_angles);
fprintf('\n========================================\n');
fprintf('  Multi-Angle SAR Imaging (%d angles)\n', N_angles);
fprintf('========================================\n');
fprintf('  Center angles: %s deg\n', num2str(center_angles));
fprintf('  Aperture:      %.1f deg\n', aperture_deg);
fprintf('  CLEAN:         maxIter=%d, thresh=%.2f, gain=%.2f\n', ...
    cleanMaxIter, cleanThreshold, cleanGain);

% 存储所有角度的散射中心（3D坐标 + 幅度 + 来源角度）
all_scatterers = [];  % columns: [x, y, z, amplitude, theta_src]

for i_ang = 1:N_angles
    theta0_deg = center_angles(i_ang);
    half_span = aperture_deg / 2;
    theta_min = theta0_deg - half_span;
    theta_max = theta0_deg + half_span;

    fprintf('\n--- Angle %d/%d: θ₀ = %.0f° [%.1f°, %.1f°] ---\n', ...
        i_ang, N_angles, theta0_deg, theta_min, theta_max);

    % --- 2a. 提取子孔径 ---
    angle_indices = find(theta_vec >= theta_min & theta_vec <= theta_max);
    if isempty(angle_indices)
        [~, sort_idx] = sort(abs(theta_vec - theta0_deg));
        n_take = max(2, floor(length(sort_idx) * aperture_deg / ...
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

    delta_theta_rad_actual = deg2rad(max(theta_sub) - min(theta_sub));
    delta_x = lambda_c / (2 * delta_theta_rad_actual);

    fprintf('  N_θ=%d, Δθ=%.4f rad (%.1f° actual, %.1f° nominal), Δx=%.3f m\n', ...
        N_theta_sub, delta_theta_rad_actual, rad2deg(delta_theta_rad_actual), aperture_deg, delta_x);

    % --- 2b. 2D 加窗 ---
    [W_2d, ~, ~] = create2DWindowSAR(N_f, N_theta_sub, windowType);
    S_win = S_sub .* W_2d;

    % --- 2c. 零填充与 2D IFFT ---
    N_r_pad = N_f * zeroPadFactor;
    N_x_pad = N_theta_sub * zeroPadFactor;
    if mod(N_r_pad, 2) ~= 0, N_r_pad = N_r_pad + 1; end
    if mod(N_x_pad, 2) ~= 0, N_x_pad = N_x_pad + 1; end

    S_padded = zeros(N_x_pad, N_r_pad);
    S_padded(1:N_theta_sub, 1:N_f) = S_win;
    I_fft = fftshift(ifft2(S_padded));

    % --- 2d. 坐标轴 ---
    delta_r_pad = delta_r * (N_f / N_r_pad);
    range_axis = (-floor(N_r_pad/2) : ceil(N_r_pad/2)-1) * delta_r_pad;
    % Cross-range 轴: 量程与孔径 Δθ 成比例。
    % IFFT 输出的 cross-range 座标系中，N_θ 个采样点对应 λ/(2·dθ) 的不模糊窗口。
    % 为避免量程受 dθ 支配而不随孔径变化，显式用孔径 Δθ 重新标定:
    %   cross_range 边界 = ± tan(Δθ/2) · (R_max/2)
    % 其中 R_max/2 取 range 轴正半轴中点，使 cross-range 与物理尺度对齐。
    R_half = max(range_axis);                        % range 轴正半轴 (m)
    X_half = R_half * tan(delta_theta_rad_actual / 2); % cross-range 半量程 (m)
    delta_x_pad = 2 * X_half / N_x_pad;              % bin 间距 (m)
    crossrange_axis = (-floor(N_x_pad/2) : ceil(N_x_pad/2)-1) * delta_x_pad;

    % --- 2e. PSF 函数 ---
    psf_func = @(r, x) sincfSAR((2*B/c) * r) .* sincfSAR((2*delta_theta_rad_actual/lambda_c) * x);

    % --- 2f. CLEAN 提取（在图像上找峰位置）---
    % 保存原始图像：CLEAN 会原地修改 I_fft 为残差
    I_orig = I_fft;
    [components, ~] = cleanAlgorithmSAR(I_fft, range_axis, crossrange_axis, ...
        psf_func, cleanThreshold, cleanMaxIter, cleanGain);

    if isempty(components)
        fprintf('  No scattering centers extracted.\n');
        continue;
    end

    % --- 2g. 用未归一化原始矩阵查表取幅度（跨角度可比）---
    % 取消 MATLAB ifft2 的 1/(Nx·Nr) 归一化，使幅度跨角度可比
    I_unnorm = abs(I_orig) * (N_x_pad * N_r_pad);
    for k = 1:size(components, 1)
        [~, r_idx] = min(abs(range_axis - components(k, 1)));
        [~, x_idx] = min(abs(crossrange_axis - components(k, 2)));
        components(k, 3) = I_unnorm(x_idx, r_idx);
    end

    fprintf('  Extracted %d scattering centers\n', size(components, 1));

    % --- 2h. 2D → 3D 坐标变换 ---
    %
    % 关键修正：散射仿真使用 exp(+j·2k·R·P) 相位约定（标准雷达 convention
    % 为 exp(-j·2k·R·P)）。MATLAB IFFT 处理该数据后，(r, x) 坐标发生
    % 符号翻转：r_sar = -k_hat·P,  x_sar = -k_perp·P。
    %
    % 因此正确重建公式为：P_3d = -r_sar·k_hat - x_sar·k_perp
    % （相位参考点在原点，重建结果直接位于 STL 模型坐标系中）
    %
    % 雷达视线方向（单站，phi=0）
    theta0_rad = deg2rad(theta0_deg);
    k_hat      = [sin(theta0_rad); 0; cos(theta0_rad)];      % Range 方向
    k_perp     = [cos(theta0_rad); 0; -sin(theta0_rad)];     % Cross-Range 方向

    % 对每个散射中心进行逆投影
    r_vals = components(:, 1);   % range (m) — IFFT 输出
    x_vals = components(:, 2);   % cross-range (m) — IFFT 输出
    A_vals = components(:, 3);   % amplitude（来自未归一化矩阵）

    % 修正后的重建公式（取反 r 和 x 以抵消相位约定引入的符号翻转）
    P_3d = - r_vals * k_hat' - x_vals * k_perp';   % N_sc × 3

    % 累积
    N_sc = size(components, 1);
    all_scatterers = [all_scatterers; ...
        P_3d, A_vals, theta0_deg * ones(N_sc, 1)];  %#ok<AGROW>

    % --- 2h. 绘制单角度 SAR 距离-方位图 (SNR 雷达图) ---
    % I_fft_dB = 20 * log10(abs(I_fft) / max(abs(I_fft(:))));
    % set(0, 'DefaultFigureVisible', 'off');
    % sarFig = figure(100 + i_ang);  % 图号偏移避免与1-3冲突
    % set(0, 'DefaultFigureVisible', 'on');
    % clf;
    % set(sarFig, 'Name', sprintf('SAR Image θ₀=%.0f°', theta0_deg), ...
    %             'NumberTitle', 'off', 'Position', [100, 100, 700, 550]);
    % imagesc(range_axis, crossrange_axis, I_fft_dB);
    % colormap('jet'); axis xy;
    % xlabel('Range (m)', 'FontSize', 11);
    % ylabel('Cross-Range (m)', 'FontSize', 11);
    % title(sprintf('SAR Image θ₀=%.0f° Δθ=%.1f° (%d CLEAN centers)', ...
    %     theta0_deg, rad2deg(delta_theta_rad_actual), N_sc), 'FontSize', 12, 'FontWeight', 'bold');
    % clim([prctile(I_fft_dB(:), 2), max(I_fft_dB(:))]);
    % colorbar;
    % hold on;
    % 标注 CLEAN 散射中心位置
    for k = 1:N_sc
        plot(components(k,1), components(k,2), 'wo', ...
            'MarkerSize', 6 + 14*components(k,3)/max(components(:,3)), ...
            'LineWidth', 1);
    end
    hold off;
end

if isempty(all_scatterers)
    error('No scattering centers extracted at any angle.');
end

fprintf('\nTotal scattering centers collected: %d\n', size(all_scatterers, 1));

%% ========================================================================
%% 3. 三维可视化
%% ========================================================================
fprintf('\nGenerating 3D visualization...\n');

% --- 参数 ---
scat_x = all_scatterers(:, 1);
scat_y = all_scatterers(:, 2);
scat_z = all_scatterers(:, 3);
scat_A = all_scatterers(:, 4);
scat_theta = all_scatterers(:, 5);

%% --- 诊断：坐标范围对比 ---
fprintf('\n--- Diagnostic: Coordinate Range Comparison ---\n');
if ~isempty(vertices)
    fprintf('  STL model bbox:\n');
    fprintf('    X: [%.3f, %.3f] m\n', min(vertices(:,1)), max(vertices(:,1)));
    fprintf('    Y: [%.3f, %.3f] m\n', min(vertices(:,2)), max(vertices(:,2)));
    fprintf('    Z: [%.3f, %.3f] m\n', min(vertices(:,3)), max(vertices(:,3)));
end
fprintf('  Scattering centers (all angles):\n');
fprintf('    X: [%.3f, %.3f] m\n', min(scat_x), max(scat_x));
fprintf('    Y: [%.3f, %.3f] m\n', min(scat_y), max(scat_y));
fprintf('    Z: [%.3f, %.3f] m\n', min(scat_z), max(scat_z));
fprintf('  Reference point P_ref: [%.3f, %.3f, %.3f] m\n', P_ref);

% 幅度转 dB，球大小与 dB 成正比（最强=0dB→200, 最弱→20）
A_dB = 20 * log10(scat_A / max(scat_A));
A_norm_dB = (scat_A - min(scat_A)) ./ max(max(scat_A - min(scat_A), 1e-10));
A_norm_4 = 1-(A_norm_dB-1).^6;
marker_sizes = 20 + 180 * A_norm_4;

% 按来源角度分配颜色
unique_thetas = unique(scat_theta);
N_colors = length(unique_thetas);
cmap = lines(N_colors);
scat_color = zeros(size(scat_theta, 1), 3);
for i = 1:length(scat_theta)
    [~, idx] = min(abs(unique_thetas - scat_theta(i)));
    scat_color(i, :) = cmap(idx, :);
end

%% --- 图1: 3D 总览（目标模型 + 散射中心）---
figure(1);
clf;
set(gcf, 'Name', '3D Scattering Centers + Target Model', ...
         'NumberTitle', 'off', 'Position', [50, 80, 1000, 750]);

hold on;

% 目标 STL 模型（半透明灰色）
if ~isempty(vertices) && ~isempty(faces)
    patch('Faces', faces, 'Vertices', vertices, ...
          'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', ...
          'FaceAlpha', 0.25, 'AmbientStrength', 0.3, ...
          'DiffuseStrength', 0.3, 'SpecularStrength', 0.1);
end

% 散射中心（颜色 = 来源角度，大小 = 幅度）
for i = 1:size(all_scatterers, 1)
    scatter3(scat_x(i), scat_y(i), scat_z(i), ...
             marker_sizes(i), scat_color(i, :), ...
             'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end

% 标注相位参考点（包围盒中心）
scatter3(P_ref(1), P_ref(2), P_ref(3), 150, 'k', 'p', ...
         'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);

% 标注各角度的雷达视线方向
scale = max(range(scat_x), max(range(scat_z))) * 0.6;
for i = 1:N_colors
    th = deg2rad(unique_thetas(i));
    k_hat = [sin(th); 0; cos(th)];
    quiver3(P_ref(1), P_ref(2), P_ref(3), ...
            scale * k_hat(1), 0, scale * k_hat(3), ...
            0, 'Color', cmap(i, :), 'LineWidth', 1.5, 'LineStyle', '--', ...
            'MaxHeadSize', 0.5, 'AutoScale', 'off');
end

hold off;
axis equal;
xlim([-1,1]); ylim([-1,1]); zlim([-1.5,3])
grid on;
box on;
xlabel('X (m)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y (m)', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Z (m)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf(['3D Scattering Center Reconstruction\n', ...
    'Model: %s | %d Angles | %d Total Centers'], ...
    modelName, N_colors, size(all_scatterers, 1)), ...
    'FontSize', 13, 'FontWeight', 'bold');

% 颜色图例
colormap(gca, cmap);
clim([min(unique_thetas), max(unique_thetas)]);
cb = colorbar;
cb.Label.String = 'Source Angle θ₀ (deg)';
cb.Label.FontSize = 11;
set(cb, 'Ticks', unique_thetas);

view(45, 25);
camlight('headlight');
lighting gouraud;
material dull;

%% --- 图2: 三视图（正视图/侧视图/俯视图）---
figure(2);
clf;
set(gcf, 'Name', 'Three-View Projections (三视图)', ...
         'NumberTitle', 'off', 'Position', [100, 100, 1100, 400]);

% 投影平面视图
view_names = {'XZ Plane (Front View)', 'YZ Plane (Side View)', 'XY Plane (Top View)'};
view_angles = {[0, 0], [90, 0], [0, 90]};
view_xlabels = {'X (m)', 'Y (m)', 'X (m)'};
view_ylabels = {'Z (m)', 'Z (m)', 'Y (m)'};
view_xdata = {scat_x, scat_y, scat_x};
view_ydata = {scat_z, scat_z, scat_y};

for v = 1:3
    subplot(1, 3, v);
    hold on;

    % 目标轮廓（投影）
    if ~isempty(vertices)
        switch v
            case 1  % XZ: front
                vx = vertices(:, 1); vy_proj = vertices(:, 3);
            case 2  % YZ: side
                vx = vertices(:, 2); vy_proj = vertices(:, 3);
            case 3  % XY: top
                vx = vertices(:, 1); vy_proj = vertices(:, 2);
        end
        % 用边界线表示目标
        k_idx = boundary(vx, vy_proj, 0.5);
        if ~isempty(k_idx)
            fill(vx(k_idx), vy_proj(k_idx), [0.85 0.85 0.85], ...
                 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 1);
        end
    end

    % 散射中心
    for i = 1:size(all_scatterers, 1)
        scatter(view_xdata{v}(i), view_ydata{v}(i), ...
                marker_sizes(i) * 0.5, scat_color(i, :), ...
                'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    end

    % 参考点
    switch v
        case 1, px = P_ref(1); py_proj = P_ref(3);
        case 2, px = P_ref(2); py_proj = P_ref(3);
        case 3, px = P_ref(1); py_proj = P_ref(2);
    end
    scatter(px, py_proj, 100, 'k', 'p', 'filled', ...
            'MarkerEdgeColor', 'w', 'LineWidth', 1);

    hold off;
    axis equal;
    grid on;
    xlabel(view_xlabels{v}, 'FontSize', 11);
    ylabel(view_ylabels{v}, 'FontSize', 11);
    title(view_names{v}, 'FontSize', 12, 'FontWeight', 'bold');
end

sgtitle(sprintf('Three-View Projections | Model: %s', modelName), ...
       'FontSize', 13, 'FontWeight', 'bold');

%% --- 图3: 散射中心幅度 vs 角度分布 ---
figure(3);
clf;
set(gcf, 'Name', 'Scattering Center Distribution (散射中心分布)', ...
         'NumberTitle', 'off', 'Position', [150, 120, 800, 350]);

subplot(1, 2, 1);
for i = 1:N_colors
    mask = (scat_theta == unique_thetas(i));
    scatter(scat_theta(mask), scat_A(mask), ...
            40 + 100*A_norm_dB(mask), cmap(i, :), ...
            'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    hold on;
end
hold off;
xlabel('Source Angle θ₀ (deg)', 'FontSize', 11);
ylabel('Scattering Amplitude', 'FontSize', 11);
title('Amplitude vs Angle', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

subplot(1, 2, 2);
histogram(scat_theta, unique_thetas, 'FaceColor', [0.3 0.5 0.8]);
xlabel('Source Angle θ₀ (deg)', 'FontSize', 11);
ylabel('Number of Centers', 'FontSize', 11);
title('Centers per Angle', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

sgtitle(sprintf('Scattering Center Statistics | Model: %s', modelName), ...
       'FontSize', 13, 'FontWeight', 'bold');

%% ========================================================================
%% 4. 保存结果
%% ========================================================================
fprintf('\nSaving results...\n');

% Create timestamped result directory
[resultDir, nowStr] = createResultDir('reconstruct_3d_scatterers');

% 保存各图
% 保存图1-3 (3D重建图)
for figNum = 1:3
    if ishandle(figNum)
        figFile = fullfile(resultDir, sprintf('scatter3d_fig%d_%s.png', figNum, nowStr));
        saveas(figNum, figFile);
        fprintf('  Figure %d: %s\n', figNum, figFile);
    end
end
% 保存各角度 SAR 距离-方位图 (figNum = 100 + i_ang)
for i_ang = 1:N_angles
    sarFig = 100 + i_ang;
    if ishandle(sarFig)
        figFile = fullfile(resultDir, sprintf('scatter3d_sar_theta%.0f_%s.png', ...
            center_angles(i_ang), nowStr));
        print(sarFig, figFile, '-dpng');
        fprintf('  SAR θ=%.0f°: %s\n', center_angles(i_ang), figFile);
    end
end

% 保存 MAT 数据
matFile = fullfile(resultDir, ['scatter3d_data_' nowStr '.mat']);
save(matFile, 'all_scatterers', 'P_ref', 'center_angles', 'unique_thetas', ...
    'lambda_c', 'B', 'c', 'delta_r', 'aperture_deg', 'modelName', '-v7.3');
fprintf('  Data: %s\n', matFile);

% 导出文本
txtFile = fullfile(resultDir, ['scatter3d_centers_' nowStr '.txt']);
fid = fopen(txtFile, 'w');
fprintf(fid, '# 3D Scattering Centers Reconstruction\n');
fprintf(fid, '# Model: %s\n', modelName);
fprintf(fid, '# Reference point: [%.4f, %.4f, %.4f] m\n', P_ref);
fprintf(fid, '# Center angles (deg): %s\n', num2str(center_angles));
fprintf(fid, '# Aperture: %.1f deg\n', aperture_deg);
fprintf(fid, '# Total centers: %d\n', size(all_scatterers, 1));
fprintf(fid, '#\n');
fprintf(fid, '# Columns: x(m)  y(m)  z(m)  amplitude  source_theta(deg)\n');
for i = 1:size(all_scatterers, 1)
    fprintf(fid, '  %.6f  %.6f  %.6f  %.6e  %.1f\n', all_scatterers(i, :));
end
fclose(fid);
fprintf('  Text: %s\n', txtFile);

fprintf('\n========================================\n');
fprintf('  3D Reconstruction Complete\n');
fprintf('========================================\n');
fprintf('  Model:            %s\n', modelName);
fprintf('  Center angles:    %s deg\n', num2str(center_angles));
fprintf('  Total centers:    %d\n', size(all_scatterers, 1));
fprintf('  Reference point:  [%.3f, %.3f, %.3f] m\n', P_ref);
fprintf('========================================\n');

end  % reconstruct_3d_scatterers

%% ========================================================================
%% 辅助函数
%% ========================================================================

% -------------------------------------------------------------------------
function [vertices, faces] = readSTLGeometry(modelFile, ~)
% 读取目标 STL 几何模型（自包含，无对话框，无副作用）
%
% 直接解析二进制 STL 文件，不依赖 stlread / stlConverter / extractCoordinatesData，
% 避免意外弹出文件选择对话框。

    vertices = [];
    faces = [];

    if isempty(modelFile) || ~exist(modelFile, 'file')
        fprintf('  STL model not found: %s (skipping model rendering)\n', modelFile);
        return;
    end

    try
        % 直接读取二进制 STL
        fid = fopen(modelFile, 'rb');
        if fid == -1
            fprintf('  Cannot open STL file: %s\n', modelFile);
            return;
        end

        % 跳过 80 字节头
        fread(fid, 80, 'uint8=>uint8');

        % 三角面元数
        nTri = fread(fid, 1, 'uint32=>uint32');

        % 预分配
        allVerts = zeros(nTri * 3, 3);

        for i = 1:nTri
            fread(fid, 3, 'float32=>double');      % 法向量（跳过）
            allVerts(3*i-2, :) = fread(fid, 3, 'float32=>double')';  % v1
            allVerts(3*i-1, :) = fread(fid, 3, 'float32=>double')';  % v2
            allVerts(3*i,   :) = fread(fid, 3, 'float32=>double')';  % v3
            fread(fid, 1, 'uint16=>uint16');       % 属性字节（跳过）
        end
        fclose(fid);

        % 去重得到唯一点集
        [vertices, ~, ic] = unique(allVerts, 'rows', 'stable');
        faces = reshape(ic, 3, nTri)';

        fprintf('  STL loaded directly: %d vertices, %d faces\n', ...
                size(vertices, 1), size(faces, 1));
    catch ME
        fprintf('  STL read error: %s\n', ME.message);
    end
end

% -------------------------------------------------------------------------
function [W_2d, w_f, w_theta] = create2DWindowSAR(N_f, N_theta, windowType)
% 二维可分离窗函数

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
function [components, residual] = cleanAlgorithmSAR(I_complex, range_axis, ...
    crossrange_axis, psf_func, threshold, maxIter, gain)
% CLEAN 迭代散射中心提取

    residual = I_complex;
    residual_mag = abs(residual);
    components = [];
    initial_peak = max(residual_mag(:));
    stop_threshold = threshold * initial_peak;

    for iter = 1:maxIter
        [peak_val, peak_idx] = max(residual_mag(:));
        if peak_val < stop_threshold, break; end

        N_r = length(range_axis);
        N_x = length(crossrange_axis);
        [peak_row, peak_col] = ind2sub([N_x, N_r], peak_idx);

        r_peak = range_axis(peak_col);
        x_peak = crossrange_axis(peak_row);
        complex_peak = residual(peak_row, peak_col);
        A_k = abs(complex_peak);
        phi_k = angle(complex_peak);

        components = [components; r_peak, x_peak, A_k, phi_k]; %#ok<AGROW>

        % 移位 PSF
        [R_mesh, X_mesh] = meshgrid(range_axis, crossrange_axis);
        psf_shifted = psf_func(R_mesh - r_peak, X_mesh - x_peak);

        residual = residual - gain * complex_peak * psf_shifted;
        residual_mag = abs(residual);
    end

    if ~isempty(components)
        fprintf('    CLEAN: %d centers, residual %.1f%%\n', ...
            size(components, 1), 100*max(residual_mag(:))/initial_peak);
    end
end

% -------------------------------------------------------------------------
function y = sincfSAR(x)
% Sinc 函数: sin(pi*x) / (pi*x)，x=0 取极限值 1

    y = ones(size(x));
    nonzero = abs(x) > 1e-15;
    y(nonzero) = sin(pi * x(nonzero)) ./ (pi * x(nonzero));
end
