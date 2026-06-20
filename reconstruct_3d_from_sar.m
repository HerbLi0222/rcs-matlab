function reconstruct_3d_from_sar(varargin)
% RECONSTRUCT_3D_FROM_SAR  基于 SAR 成像结果的三维散射中心重建
%
%   使用 main_sar_imaging.m 的 SAR 处理流程（2D IFFT + CLEAN + PSF），
%   在多个中心观测角度下分别成像，将各角度的 (range, cross-range)
%   散射中心二维坐标逆投影回三维空间，并与 STL 模型叠加显示。
%
%   与 reconstruct_3d_scatterers.m 的区别：
%     - reconstruct_3d_scatterers: 简化版 SAR 处理（单一孔径参数）
%     - reconstruct_3d_from_sar:     使用 main_sar_imaging 的完整流程
%       （多子孔径、FFT vs PSF 对比、展宽测量等），然后多角度 3D 重建
%
%   方法原理:
%     散射仿真使用 exp(+j·2k·R·P) 相位约定（标准雷达用 exp(-j·2k·R·P)），
%     MATLAB IFFT 后 (r, x) 坐标发生符号翻转：
%       r_sar = -k̂·P,   x_sar = -k̂_perp·P
%
%     修正后的逆变换（从2D到3D）：
%       P_3D = - r_sar·k̂ - x_sar·k̂_perp
%
%     其中：
%       k̂       = [sin(θ₀), 0, cos(θ₀)]   — 雷达视线方向 (Range)
%       k̂_perp  = [cos(θ₀), 0, -sin(θ₀)]  — 方位向 (Cross-Range)
%
%   Usage:
%     >> reconstruct_3d_from_sar
%     >> reconstruct_3d_from_sar('DataFile', 'wideband_scattering_xxx.mat')
%     >> reconstruct_3d_from_sar('CenterAngles', [30, 60, 90, 120, 150])
%     >> reconstruct_3d_from_sar('SubAperturesDeg', [10, 5], ...
%                                'CleanMaxIter', 50, 'CleanThreshold', 0.03)

%% ========================================================================
%% 0. 参数解析
%% ========================================================================
addpath('lib');

p = inputParser;
p.addOptional('DataFile', '', @(x) ischar(x) || isstring(x) || isempty(x));
p.addParameter('CenterAngles', 30:15:120, @(x) isnumeric(x));
p.addParameter('SubAperturesDeg', 10, @(x) isnumeric(x) && isscalar(x));
p.addParameter('WindowType', 'kaiser', @(x) ischar(x) || isstring(x));
p.addParameter('ZeroPadFactor', 2, @(x) isnumeric(x) && isscalar(x));
p.addParameter('CleanMaxIter', 50, @(x) isnumeric(x) && isscalar(x));
p.addParameter('CleanThreshold', 0.05, @(x) isnumeric(x) && isscalar(x));
p.addParameter('CleanGain', 0.5, @(x) isnumeric(x) && isscalar(x));

if ~isempty(varargin) && (ischar(varargin{1}) || isstring(varargin{1})) ...
        && any(strcmp(char(varargin{1}), {'DataFile','CenterAngles',...
        'SubAperturesDeg','WindowType','ZeroPadFactor','CleanMaxIter',...
        'CleanThreshold','CleanGain'}))
    varargin = [{''}, varargin];
end
p.parse(varargin{:});

resultFile        = p.Results.DataFile;
center_angles_deg = p.Results.CenterAngles;
subApertureDeg    = p.Results.SubAperturesDeg;
windowType        = p.Results.WindowType;
zeroPadFactor     = p.Results.ZeroPadFactor;
cleanMaxIter      = p.Results.CleanMaxIter;
cleanThreshold    = p.Results.CleanThreshold;
cleanGain         = p.Results.CleanGain;

%% ========================================================================
%% 1. 加载宽带散射数据和目标几何
%% ========================================================================
fprintf('========================================\n');
fprintf('  3D Reconstruction from SAR Imaging\n');
fprintf('  (基于 SAR 成像结果的三维重建)\n');
fprintf('========================================\n\n');

% --- 1a. 加载 wideband_scattering 数据 ---
if isempty(resultFile)
    resultFile = findLatestResultFile('wideband_scattering_*.mat');
    if isempty(resultFile)
        error('No wideband_scattering_*.mat found in results/.\nRun main_wideband_scattering first.');
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

% 模型信息
if isfield(data, 'inputModel')
    [~, modelName, ~] = fileparts(data.inputModel);
    modelFile = fullfile('stl_models', data.inputModel);
else
    modelName = 'Unknown';
    modelFile = '';
end

% 相位参考点（仅用于可视化标注，不用于重建公式）
if isfield(data, 'bbox_center')
    P_ref = data.bbox_center(:);
else
    P_ref = [0; 0; 0];
end

% 构建角度向量
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

fprintf('  Model: %s, Bandwidth: %.1f GHz, λc: %.4f m\n', ...
    modelName, B/1e9, lambda_c);
fprintf('  Freq points: %d, Angle range: [%.1f, %.1f]° (%d angles)\n', ...
    N_f, theta_vec(1), theta_vec(end), length(theta_vec));
fprintf('  Reference point (visual only): [%.3f, %.3f, %.3f] m\n', P_ref);

% --- 1b. 加载 STL 模型 ---
fprintf('\nLoading target geometry...\n');
[vertices, faces] = readSTLGeometry(modelFile);
fprintf('  Vertices: %d, Faces: %d\n', size(vertices, 1), size(faces, 1));

%% ========================================================================
%% 2. 多角度 SAR 成像（使用 main_sar_imaging 流程）
%% ========================================================================
N_angles = length(center_angles_deg);
fprintf('\n========================================\n');
fprintf('  Multi-Angle SAR Imaging (%d center angles)\n', N_angles);
fprintf('========================================\n');
fprintf('  Center angles:    %s deg\n', num2str(center_angles_deg));
fprintf('  Sub-aperture:     %.1f deg\n', subApertureDeg);
fprintf('  Window:           %s\n', windowType);
fprintf('  Zero-pad factor:  %dx\n', zeroPadFactor);
fprintf('  CLEAN:            maxIter=%d, thresh=%.2f, gain=%.2f\n', ...
    cleanMaxIter, cleanThreshold, cleanGain);

% 存储所有角度的散射中心（3D坐标 + 幅度 + 来源角度）
all_scatterers = [];  % columns: [x, y, z, amplitude, theta_src]

% 角度循环
for i_ang = 1:N_angles
    theta0_deg = center_angles_deg(i_ang);
    half_span = subApertureDeg / 2;
    theta_min = theta0_deg - half_span;
    theta_max = theta0_deg + half_span;

    fprintf('\n--- Center Angle %d/%d: θ₀ = %.1f° [%.1f°, %.1f°] ---\n', ...
        i_ang, N_angles, theta0_deg, theta_min, theta_max);

    % --- 2a. 提取子孔径数据 ---
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

    % 使用实际角度跨度进行方位分辨率校准
    delta_theta_rad_actual = deg2rad(max(theta_sub) - min(theta_sub));
    delta_x = lambda_c / (2 * delta_theta_rad_actual);

    fprintf('  N_θ=%d, Δθ=%.4f rad (%.1f° actual, %.1f° nominal), Δx=%.3f m\n', ...
        N_theta_sub, delta_theta_rad_actual, rad2deg(delta_theta_rad_actual), ...
        subApertureDeg, delta_x);

    % --- 2b. 零填充参数 ---
    N_r_pad = N_f * zeroPadFactor;
    N_x_pad = N_theta_sub * zeroPadFactor;
    if mod(N_r_pad, 2) ~= 0, N_r_pad = N_r_pad + 1; end
    if mod(N_x_pad, 2) ~= 0, N_x_pad = N_x_pad + 1; end

    % --- 2c. 二维加窗 ---
    [W_2d, ~, ~] = create2DWindow(N_f, N_theta_sub, windowType);
    S_win = S_sub .* W_2d;

    % --- 2d. 补零与 2D IFFT（FFT 法成像）---
    S_padded = padarray(S_win, [N_x_pad - N_theta_sub, N_r_pad - N_f], 0, 'post');
    I_fft = fftshift(ifft2(S_padded));

    % --- 2e. 坐标轴标定 ---
    delta_r_pad = delta_r * (N_f / N_r_pad);
    range_axis = (-floor(N_r_pad/2) : ceil(N_r_pad/2)-1) * delta_r_pad;
    % Cross-range 轴: 量程与孔径 Δθ 成比例
    R_half = max(range_axis);
    X_half = R_half * tan(delta_theta_rad_actual / 2);
    delta_x_pad = 2 * X_half / N_x_pad;
    crossrange_axis = (-floor(N_x_pad/2) : ceil(N_x_pad/2)-1) * delta_x_pad;

    % --- 2f. PSF 函数 ---
    psf_func = @(r, x) sincf((2*B/c) * r) .* sincf((2*delta_theta_rad_actual/lambda_c) * x);

    % --- 2g. CLEAN 算法提取散射中心 ---
    [clean_components, ~] = cleanAlgorithm(I_fft, range_axis, crossrange_axis, ...
        psf_func, cleanThreshold, cleanMaxIter, cleanGain);

    if isempty(clean_components)
        fprintf('  No scattering centers extracted.\n');
        continue;
    end

    N_clean = size(clean_components, 1);
    fprintf('  CLEAN extracted %d scattering centers\n', N_clean);

    % --- 2h. PSF 叠加合成图像（用于可视化）---
    [I_psf, ~] = buildPSFImage(clean_components, range_axis, crossrange_axis, psf_func);
    I_psf_dB = 20 * log10(abs(I_psf) / max(abs(I_psf(:))));
    I_fft_dB = 20 * log10(abs(I_fft) / max(abs(I_fft(:))));

    % --- 2i. 2D → 3D 坐标变换（修正符号翻转）---
    %
    % 关键修正：散射仿真使用 exp(+j·2k·R·P) 相位约定（标准雷达 convention
    % 为 exp(-j·2k·R·P)）。MATLAB IFFT 处理该数据后，(r, x) 坐标发生
    % 符号翻转：r_sar = -k_hat·P,  x_sar = -k_perp·P。
    %
    % 因此正确重建公式为：P_3d = -r_sar·k_hat - x_sar·k_perp
    %
    theta0_rad = deg2rad(theta0_deg);
    k_hat  = [sin(theta0_rad); 0; cos(theta0_rad)];      % Range 方向
    k_perp = [cos(theta0_rad); 0; -sin(theta0_rad)];     % Cross-Range 方向

    r_vals = clean_components(:, 1);   % range (m) — IFFT 输出
    x_vals = clean_components(:, 2);   % cross-range (m) — IFFT 输出
    A_vals = clean_components(:, 3);   % amplitude

    % 修正后的重建公式（取反 r 和 x）
    P_3d = - r_vals * k_hat' - x_vals * k_perp';   % N_sc × 3

    % 累积
    all_scatterers = [all_scatterers; ...
        P_3d, A_vals, theta0_deg * ones(N_clean, 1)];  %#ok<AGROW>

    % --- 2j. 单角度 SAR 图像快照（可选：关闭以加速）---
    if N_angles <= 5
        plotSingleAngleSAR(I_fft_dB, I_psf_dB, range_axis, crossrange_axis, ...
            clean_components, theta0_deg, subApertureDeg, modelName);
    end

end  % 角度循环

if isempty(all_scatterers)
    error('No scattering centers extracted at any angle.');
end

fprintf('\nTotal scattering centers collected: %d\n', size(all_scatterers, 1));

%% ========================================================================
%% 3. 三维可视化
%% ========================================================================
fprintf('\nGenerating 3D visualization...\n');

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
fprintf('  Scattering centers (all angles, corrected):\n');
fprintf('    X: [%.3f, %.3f] m\n', min(scat_x), max(scat_x));
fprintf('    Y: [%.3f, %.3f] m\n', min(scat_y), max(scat_y));
fprintf('    Z: [%.3f, %.3f] m\n', min(scat_z), max(scat_z));

% 幅度归一化（用于点大小）
A_norm = scat_A / max(scat_A);
marker_sizes = 20 + 180 * A_norm;   % 20 ~ 200

% 按来源角度着色
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
set(gcf, 'Name', '3D Scattering Centers from SAR (SAR 三维散射中心)', ...
         'NumberTitle', 'off', 'Position', [50, 80, 1000, 750]);

hold on;

% STL 模型（半透明灰色）
if ~isempty(vertices) && ~isempty(faces)
    patch('Faces', faces, 'Vertices', vertices, ...
          'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', ...
          'FaceAlpha', 0.25, 'AmbientStrength', 0.3, ...
          'DiffuseStrength', 0.3, 'SpecularStrength', 0.1);
end

% 散射中心
for i = 1:size(all_scatterers, 1)
    scatter3(scat_x(i), scat_y(i), scat_z(i), ...
             marker_sizes(i), scat_color(i, :), ...
             'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end

% 参考点（仅用于可视化标注）
scatter3(P_ref(1), P_ref(2), P_ref(3), 150, 'k', 'p', ...
         'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);

% 各角度雷达视线方向
scale = max(range(scat_x), range(scat_z)) * 0.6;
for i = 1:N_colors
    th = deg2rad(unique_thetas(i));
    k_dir = [sin(th); 0; cos(th)];
    quiver3(P_ref(1), P_ref(2), P_ref(3), ...
            scale * k_dir(1), 0, scale * k_dir(3), ...
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
title(sprintf(['3D Scattering Centers from SAR Imaging\n', ...
    'Model: %s | %d Center Angles | %d Total Centers'], ...
    modelName, N_colors, size(all_scatterers, 1)), ...
    'FontSize', 13, 'FontWeight', 'bold');

colormap(gca, cmap);
clim([min(unique_thetas), max(unique_thetas)]);
cb = colorbar;
cb.Label.String = 'Center Angle θ₀ (deg)';
cb.Label.FontSize = 11;
set(cb, 'Ticks', unique_thetas);

view(45, 25);
camlight('headlight');
lighting gouraud;
material dull;

%% --- 图2: 三视图 ---
figure(2);
clf;
set(gcf, 'Name', 'Three-View Projections (三视图)', ...
         'NumberTitle', 'off', 'Position', [100, 100, 1100, 400]);

view_configs = {
    'XZ Plane (Front)',  'X (m)', 'Z (m)', scat_x, scat_z;
    'YZ Plane (Side)',   'Y (m)', 'Z (m)', scat_y, scat_z;
    'XY Plane (Top)',    'X (m)', 'Y (m)', scat_x, scat_y;
};

for v = 1:3
    subplot(1, 3, v);
    hold on;

    % STL 轮廓投影
    if ~isempty(vertices)
        switch v
            case 1, bx = vertices(:,1); by_p = vertices(:,3);
            case 2, bx = vertices(:,2); by_p = vertices(:,3);
            case 3, bx = vertices(:,1); by_p = vertices(:,2);
        end
        try
            k_idx = boundary(bx, by_p, 0.5);
            if ~isempty(k_idx)
                fill(bx(k_idx), by_p(k_idx), [0.85 0.85 0.85], ...
                     'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 1);
            end
        catch
        end
    end

    % 散射中心
    for i = 1:size(all_scatterers, 1)
        scatter(view_configs{v,4}(i), view_configs{v,5}(i), ...
                marker_sizes(i) * 0.5, scat_color(i, :), ...
                'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    end

    % 参考点
    switch v
        case 1, px = P_ref(1); py_p = P_ref(3);
        case 2, px = P_ref(2); py_p = P_ref(3);
        case 3, px = P_ref(1); py_p = P_ref(2);
    end
    scatter(px, py_p, 100, 'k', 'p', 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1);

    hold off;
    axis equal;
    grid on;
    xlabel(view_configs{v,2}, 'FontSize', 11);
    ylabel(view_configs{v,3}, 'FontSize', 11);
    title(view_configs{v,1}, 'FontSize', 12, 'FontWeight', 'bold');
end

sgtitle(sprintf('Three-View Projections | Model: %s', modelName), ...
       'FontSize', 13, 'FontWeight', 'bold');

%% --- 图3: 散射中心幅度 vs 角度 ---
figure(3);
clf;
set(gcf, 'Name', 'Scattering Center Distribution (散射中心分布)', ...
         'NumberTitle', 'off', 'Position', [150, 120, 800, 350]);

subplot(1, 2, 1);
hold on;
for i = 1:N_colors
    mask = (scat_theta == unique_thetas(i));
    scatter(scat_theta(mask), scat_A(mask), ...
            40 + 100*A_norm(mask), cmap(i, :), ...
            'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
end
hold off;
xlabel('Center Angle θ₀ (deg)', 'FontSize', 11);
ylabel('Scattering Amplitude', 'FontSize', 11);
title('Amplitude vs Angle', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

subplot(1, 2, 2);
histogram(scat_theta, unique_thetas, 'FaceColor', [0.3 0.5 0.8]);
xlabel('Center Angle θ₀ (deg)', 'FontSize', 11);
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
[resultDir, nowStr] = createResultDir('reconstruct_3d_from_sar');

% 保存各图 (图1-3: 3D重建图)
for figNum = 1:3
    if ishandle(figNum)
        figFile = fullfile(resultDir, sprintf('sar3d_fig%d_%s.png', figNum, nowStr));
        saveas(figNum, figFile);
        fprintf('  Figure %d: %s\n', figNum, figFile);
    end
end
% 保存各角度 SAR 距离-方位图 (figNum = 100 + θ₀)
for i_ang = 1:N_angles
    sarFigNum = 100 + center_angles_deg(i_ang);
    if ishandle(sarFigNum)
        figFile = fullfile(resultDir, sprintf('sar3d_sar_theta%.0f_%s.png', ...
            center_angles_deg(i_ang), nowStr));
        print(sarFigNum, figFile, '-dpng');
        fprintf('  SAR θ=%.0f°: %s\n', center_angles_deg(i_ang), figFile);
    end
end

% 保存 MAT 数据
matFile = fullfile(resultDir, ['sar3d_data_' nowStr '.mat']);
save(matFile, 'all_scatterers', 'P_ref', 'center_angles_deg', 'unique_thetas', ...
    'subApertureDeg', 'lambda_c', 'B', 'c', 'delta_r', 'windowType', ...
    'cleanThreshold', 'cleanMaxIter', 'cleanGain', 'modelName', '-v7.3');
fprintf('  Data: %s\n', matFile);

% 导出文本
txtFile = fullfile(resultDir, ['sar3d_centers_' nowStr '.txt']);
fid = fopen(txtFile, 'w');
fprintf(fid, '# 3D Scattering Centers from SAR Imaging\n');
fprintf(fid, '# Model: %s\n', modelName);
fprintf(fid, '# Center angles (deg): %s\n', num2str(center_angles_deg));
fprintf(fid, '# Sub-aperture: %.1f deg\n', subApertureDeg);
fprintf(fid, '# Window: %s, CLEAN: maxIter=%d, thresh=%.2f, gain=%.2f\n', ...
    windowType, cleanMaxIter, cleanThreshold, cleanGain);
fprintf(fid, '# Total centers: %d\n', size(all_scatterers, 1));
fprintf(fid, '# Phase convention corrected: P_3d = -r*k_hat - x*k_perp\n');
fprintf(fid, '#\n');
fprintf(fid, '# Columns: x(m)  y(m)  z(m)  amplitude  source_theta(deg)\n');
for i = 1:size(all_scatterers, 1)
    fprintf(fid, '  %.6f  %.6f  %.6f  %.6e  %.1f\n', all_scatterers(i, :));
end
fclose(fid);
fprintf('  Text: %s\n', txtFile);

fprintf('\n========================================\n');
fprintf('  SAR-Based 3D Reconstruction Complete\n');
fprintf('========================================\n');
fprintf('  Model:              %s\n', modelName);
fprintf('  Center angles:      %s deg\n', num2str(center_angles_deg));
fprintf('  Sub-aperture:       %.1f deg\n', subApertureDeg);
fprintf('  Total centers:      %d\n', size(all_scatterers, 1));
fprintf('  Phase correction:   P_3d = -r·k̂ - x·k̂_perp\n');
fprintf('========================================\n');

end  % reconstruct_3d_from_sar

%% ========================================================================
%% 辅助函数
%% ========================================================================

% -------------------------------------------------------------------------
function [vertices, faces] = readSTLGeometry(modelFile)
% 直接读取二进制 STL 文件

    vertices = [];
    faces = [];

    if isempty(modelFile) || ~exist(modelFile, 'file')
        fprintf('  STL model not found: %s (skipping model rendering)\n', modelFile);
        return;
    end

    try
        fid = fopen(modelFile, 'rb');
        if fid == -1
            fprintf('  Cannot open STL file: %s\n', modelFile);
            return;
        end

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

        fprintf('  STL loaded directly: %d vertices, %d faces\n', ...
                size(vertices, 1), size(faces, 1));
    catch ME
        fprintf('  STL read error: %s\n', ME.message);
    end
end

% -------------------------------------------------------------------------
function [W_2d, w_f, w_theta] = create2DWindow(N_f, N_theta, windowType)
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
function [I_psf, I_psf_mag] = buildPSFImage(components, range_axis, ...
    crossrange_axis, psf_func)
% 基于 CLEAN 提取的散射中心构建 PSF 叠加图像

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
function plotSingleAngleSAR(I_fft_dB, I_psf_dB, range_axis, crossrange_axis, ...
    components, theta0_deg, subApertureDeg, modelName)
% 为单个中心角度绘制 FFT vs PSF 对比图

    N_clean = size(components, 1);
    figNum = 100 + theta0_deg;  % 为每个角度分配唯一图号

    set(0, 'DefaultFigureVisible', 'off');
    figure(figNum);
    set(0, 'DefaultFigureVisible', 'on');
    clf;
    set(gcf, 'Name', sprintf('SAR Image θ₀=%.0f°', theta0_deg), ...
             'NumberTitle', 'off', 'Position', [200, 100, 900, 450]);

    % 全局色彩范围
    all_vals = [I_fft_dB(:); I_psf_dB(:)];
    c_min = prctile(all_vals, 1);
    c_max = max(all_vals);

    % FFT 图像
    subplot(1, 2, 1);
    imagesc(range_axis, crossrange_axis, I_fft_dB);
    colormap('jet');
    axis xy;
    xlabel('Range (m)'); ylabel('Cross-Range (m)');
    title(sprintf('FFT SAR θ₀=%.0f° Δθ=%.1f°', theta0_deg, subApertureDeg));
    caxis([c_min, c_max]); colorbar;

    % PSF + CLEAN 图像
    subplot(1, 2, 2);
    hold on;
    imagesc(range_axis, crossrange_axis, I_psf_dB);
    colormap('jet');
    axis xy;
    if ~isempty(components)
        amps = components(:, 3);
        amp_norm = amps / max(amps);
        for k = 1:N_clean
            plot(components(k, 1), components(k, 2), 'wo', ...
                'MarkerSize', 5 + 10*amp_norm(k), 'LineWidth', 1.2);
        end
    end
    xlabel('Range (m)'); ylabel('Cross-Range (m)');
    title(sprintf('PSF+CLEAN θ₀=%.0f° (%d centers)', theta0_deg, N_clean));
    caxis([c_min, c_max]); colorbar;
    hold off;

    sgtitle(sprintf('SAR Imaging | Model: %s | θ₀=%.0f°', modelName, theta0_deg), ...
           'FontSize', 12, 'FontWeight', 'bold');
end

% -------------------------------------------------------------------------
function y = sincf(x)
% Sinc 函数: sin(pi*x) / (pi*x)

    y = ones(size(x));
    nonzero = abs(x) > 1e-15;
    y(nonzero) = sin(pi * x(nonzero)) ./ (pi * x(nonzero));
end
