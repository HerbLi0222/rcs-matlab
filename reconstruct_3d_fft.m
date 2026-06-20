function reconstruct_3d_fft(varargin)
% RECONSTRUCT_3D_FFT  FFT 峰值检测的三维散射中心重建
%
%   从 wideband_scattering 数据中，在多个中心观测角度下进行 SAR 成像，
%   使用 FFT 法（2D IFFT + 峰值检测）提取散射中心，逆投影回三维空间，
%   并与目标 STL 模型叠加显示。
%
%   与 reconstruct_3d_scatterers.m 的区别:
%     - reconstruct_3d_scatterers: CLEAN + PSF 叠加
%     - reconstruct_3d_fft:       FFT 峰值检测（无 CLEAN 迭代）
%
%   Usage:
%     >> reconstruct_3d_fft
%     >> reconstruct_3d_fft('DataFile', 'wideband_scattering_xxx.mat')
%     >> reconstruct_3d_fft('CenterAngles', [30, 60, 90, 120, 150], ...
%           'ApertureDeg', 10, 'FftThreshold', -20)

%% ========================================================================
%% 0. 参数解析
%% ========================================================================
addpath('lib');

p = inputParser;
p.addOptional('DataFile', '', @(x) ischar(x) || isstring(x) || isempty(x));
p.addParameter('CenterAngles', 0:15:360, @(x) isnumeric(x));
p.addParameter('ApertureDeg', 10, @(x) isnumeric(x) && isscalar(x));
p.addParameter('WindowType', 'chebwin', @(x) ischar(x) || isstring(x));
p.addParameter('ZeroPadFactor', 2, @(x) isnumeric(x) && isscalar(x));
p.addParameter('FftThreshold', -20, @(x) isnumeric(x) && isscalar(x));   % dB 阈值
p.addParameter('FftMinSep', 0.04, @(x) isnumeric(x) && isscalar(x));     % 最小峰间距 (m)
p.addParameter('FftMaxPeaks', 100, @(x) isnumeric(x) && isscalar(x));    % 每角度最大峰数

if ~isempty(varargin) && (ischar(varargin{1}) || isstring(varargin{1})) ...
        && any(strcmp(char(varargin{1}), {'DataFile','CenterAngles',...
        'ApertureDeg','WindowType','ZeroPadFactor','FftThreshold',...
        'FftMinSep','FftMaxPeaks'}))
    varargin = [{''}, varargin];
end
p.parse(varargin{:});

resultFile    = p.Results.DataFile;
center_angles = p.Results.CenterAngles;
aperture_deg  = p.Results.ApertureDeg;
windowType    = p.Results.WindowType;
zeroPadFactor = p.Results.ZeroPadFactor;
fftThreshold  = p.Results.FftThreshold;
fftMinSep     = p.Results.FftMinSep;
fftMaxPeaks   = p.Results.FftMaxPeaks;

%% ========================================================================
%% 1. 加载宽带散射数据和目标几何
%% ========================================================================
fprintf('========================================\n');
fprintf('  FFT-Based 3D Scattering Center Reconstruction\n');
fprintf('  (基于 FFT 峰值检测的三维散射中心重建)\n');
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
fprintf('  Freq points: %d, Angle range: [%.0f, %.0f]°\n', N_f, theta_vec(1), theta_vec(end));

% --- 1b. 加载 STL 模型 ---
fprintf('\nLoading target geometry...\n');
[vertices, faces] = readSTLGeometry(modelFile, data);
fprintf('  Vertices: %d, Faces: %d\n', size(vertices, 1), size(faces, 1));

%% ========================================================================
%% 2. 多角度 FFT SAR 成像与峰值检测
%% ========================================================================
N_angles = length(center_angles);
fprintf('\n========================================\n');
fprintf('  Multi-Angle FFT SAR Imaging (%d angles)\n', N_angles);
fprintf('========================================\n');
fprintf('  Center angles:  %s deg\n', num2str(center_angles));
fprintf('  Aperture:       %.1f deg\n', aperture_deg);
fprintf('  FFT threshold:  %.0f dB\n', fftThreshold);
fprintf('  Min separation: %.2f m\n', fftMinSep);

all_scatterers = [];  % [x, y, z, amplitude, theta_src]

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
    R_half = max(range_axis);
    X_half = R_half * tan(delta_theta_rad_actual / 2);
    delta_x_pad = 2 * X_half / N_x_pad;
    crossrange_axis = (-floor(N_x_pad/2) : ceil(N_x_pad/2)-1) * delta_x_pad;

    % --- 2f. FFT 峰值检测 ---
    I_fft_mag = abs(I_fft);
    I_fft_dB = 20 * log10(I_fft_mag / max(I_fft_mag(:)));

    [fft_peaks, ~] = detectPeaksFFT(I_fft_mag, I_fft_dB, range_axis, ...
        crossrange_axis, fftThreshold, fftMinSep, fftMaxPeaks);

    N_peaks = size(fft_peaks, 1);
    fprintf('  FFT peaks detected: %d\n', N_peaks);

    % --- 2g. 2D → 3D 逆投影 ---
    theta0_rad = deg2rad(theta0_deg);
    k_hat  = [sin(theta0_rad); 0; cos(theta0_rad)];
    k_perp = [cos(theta0_rad); 0; -sin(theta0_rad)];

    if ~isempty(fft_peaks)
        r_f = fft_peaks(:, 1);
        x_f = fft_peaks(:, 2);
        A_f = fft_peaks(:, 3);
        P_3d = - r_f * k_hat' - x_f * k_perp';
        all_scatterers = [all_scatterers; P_3d, A_f, theta0_deg * ones(N_peaks, 1)]; %#ok<AGROW>
    end

    % --- 2h. 单角度 SAR 距离-方位图 ---
    % set(0, 'DefaultFigureVisible', 'off');
    % sarFig = figure(100 + i_ang);
    % set(0, 'DefaultFigureVisible', 'on');
    % clf;
    % set(sarFig, 'Name', sprintf('FFT SAR Image θ₀=%.0f°', theta0_deg), ...
    %             'NumberTitle', 'off', 'Position', [100, 100, 700, 550]);
    % imagesc(range_axis, crossrange_axis, I_fft_dB);
    % colormap('jet'); axis xy;
    % xlabel('Range (m)', 'FontSize', 11);
    % ylabel('Cross-Range (m)', 'FontSize', 11);
    % title(sprintf('FFT SAR θ₀=%.0f° Δθ=%.1f° (%d peaks)', ...
    %     theta0_deg, rad2deg(delta_theta_rad_actual), N_peaks), 'FontSize', 12, 'FontWeight', 'bold');
    % caxis([prctile(I_fft_dB(:), 2), max(I_fft_dB(:))]);
    % colorbar;
    % hold on;
    % if ~isempty(fft_peaks)
    %     peak_sizes = 20 + 80 * fft_peaks(:,3) / max(fft_peaks(:,3));
    %     scatter(fft_peaks(:,1), fft_peaks(:,2), peak_sizes, 'wo', ...
    %         'LineWidth', 1);
    % end
    % hold off;
end

if isempty(all_scatterers)
    error('No scattering centers extracted at any angle.');
end

fprintf('\nTotal FFT scattering centers: %d\n', size(all_scatterers, 1));

%% ========================================================================
%% 3. 三维可视化
%% ========================================================================
fprintf('\nGenerating 3D visualization...\n');

scat_x = all_scatterers(:, 1);
scat_y = all_scatterers(:, 2);
scat_z = all_scatterers(:, 3);
scat_A = all_scatterers(:, 4);
scat_theta = all_scatterers(:, 5);

% 诊断：坐标范围
fprintf('\n--- Coordinate Range Comparison ---\n');
if ~isempty(vertices)
    fprintf('  STL bbox:\n');
    fprintf('    X: [%.3f, %.3f]\n', min(vertices(:,1)), max(vertices(:,1)));
    fprintf('    Y: [%.3f, %.3f]\n', min(vertices(:,2)), max(vertices(:,2)));
    fprintf('    Z: [%.3f, %.3f]\n', min(vertices(:,3)), max(vertices(:,3)));
end
fprintf('  FFT scatterers:\n');
fprintf('    X: [%.3f, %.3f]\n', min(scat_x), max(scat_x));
fprintf('    Y: [%.3f, %.3f]\n', min(scat_y), max(scat_y));
fprintf('    Z: [%.3f, %.3f]\n', min(scat_z), max(scat_z));

% 幅度归一化
A_norm = scat_A / max(scat_A);
marker_sizes = 40 + 160 * A_norm;
A_norm_dB = (scat_A - min(scat_A)) ./ max(max(scat_A - min(scat_A), 1e-10));
A_norm_4 = 1-(A_norm_dB-1).^8;
marker_sizes = 20 + 180 * A_norm_4;

% 按角度分配颜色
unique_thetas = unique(scat_theta);
N_colors = length(unique_thetas);
cmap = lines(N_colors);
scat_color = zeros(length(scat_theta), 3);
for i = 1:length(scat_theta)
    [~, idx] = min(abs(unique_thetas - scat_theta(i)));
    scat_color(i, :) = cmap(idx, :);
end

%% --- 图1: 3D 总览 ---
figure(1);
clf;
set(gcf, 'Name', 'FFT 3D Scattering Centers + Target Model', ...
         'NumberTitle', 'off', 'Position', [50, 80, 1000, 750]);
hold on;

if ~isempty(vertices) && ~isempty(faces)
    patch('Faces', faces, 'Vertices', vertices, ...
          'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', ...
          'FaceAlpha', 0.25, 'AmbientStrength', 0.3, ...
          'DiffuseStrength', 0.3, 'SpecularStrength', 0.1);
end

for i = 1:size(all_scatterers, 1)
    scatter3(scat_x(i), scat_y(i), scat_z(i), ...
             marker_sizes(i), scat_color(i, :), ...
             'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end

scatter3(P_ref(1), P_ref(2), P_ref(3), 150, 'k', 'p', 'filled', ...
         'MarkerEdgeColor', 'w', 'LineWidth', 1.5);

% 标注各角度雷达视线方向
scale = max(range(scat_x), range(scat_z)) * 0.6;
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
grid on; box on;
xlim([-1,1]), ylim([-1,1]), zlim([-1.5, 3])
xlabel('X (m)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y (m)', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Z (m)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('FFT 3D Scattering Centers | Model: %s | %d Angles | %d Centers', ...
    modelName, N_colors, size(all_scatterers, 1)), 'FontSize', 13, 'FontWeight', 'bold');

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

%% --- 图2: 三视图 ---
figure(2);
clf;
set(gcf, 'Name', 'FFT Three-View Projections', ...
         'NumberTitle', 'off', 'Position', [100, 100, 1100, 400]);

view_names = {'XZ Plane (Front)', 'YZ Plane (Side)', 'XY Plane (Top)'};
view_angles = {[0, 0], [90, 0], [0, 90]};
view_xlabels = {'X (m)', 'Y (m)', 'X (m)'};
view_ylabels = {'Z (m)', 'Z (m)', 'Y (m)'};
view_xdata = {scat_x, scat_y, scat_x};
view_ydata = {scat_z, scat_z, scat_y};

for v = 1:3
    subplot(1, 3, v);
    hold on;

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

    for i = 1:size(all_scatterers, 1)
        scatter(view_xdata{v}(i), view_ydata{v}(i), ...
                marker_sizes(i) * 0.5, scat_color(i, :), ...
                'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    end

    switch v
        case 1, px = P_ref(1); py_p = P_ref(3);
        case 2, px = P_ref(2); py_p = P_ref(3);
        case 3, px = P_ref(1); py_p = P_ref(2);
    end
    scatter(px, py_p, 100, 'k', 'p', 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1);

    hold off;
    axis equal;
    grid on;
    xlabel(view_xlabels{v}, 'FontSize', 11);
    ylabel(view_ylabels{v}, 'FontSize', 11);
    title(view_names{v}, 'FontSize', 12, 'FontWeight', 'bold');
end

sgtitle(sprintf('FFT Three-View Projections | Model: %s', modelName), ...
       'FontSize', 13, 'FontWeight', 'bold');

%% --- 图3: 散射中心幅度 vs 角度 ---
figure(3);
clf;
set(gcf, 'Name', 'FFT Scattering Center Distribution', ...
         'NumberTitle', 'off', 'Position', [150, 120, 800, 350]);

subplot(1, 2, 1);
for i = 1:N_colors
    mask = (scat_theta == unique_thetas(i));
    scatter(scat_theta(mask), scat_A(mask), ...
            40 + 100*A_norm(mask), cmap(i, :), ...
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

sgtitle(sprintf('FFT Scattering Center Statistics | Model: %s', modelName), ...
       'FontSize', 13, 'FontWeight', 'bold');

%% ========================================================================
%% 4. 保存结果
%% ========================================================================
fprintf('\nSaving results...\n');

[resultDir, nowStr] = createResultDir('reconstruct_3d_fft');

% 图1-3
for figNum = 1:3
    if ishandle(figNum)
        figFile = fullfile(resultDir, sprintf('fft3d_fig%d_%s.png', figNum, nowStr));
        saveas(figNum, figFile);
        fprintf('  Figure %d: %s\n', figNum, figFile);
    end
end

% 各角度 SAR 图
for i_ang = 1:N_angles
    sarFig = 100 + i_ang;
    if ishandle(sarFig)
        figFile = fullfile(resultDir, sprintf('fft3d_sar_theta%.0f_%s.png', ...
            center_angles(i_ang), nowStr));
        print(sarFig, figFile, '-dpng');
        fprintf('  SAR θ=%.0f°: %s\n', center_angles(i_ang), figFile);
    end
end

% MAT 数据
matFile = fullfile(resultDir, ['fft3d_data_' nowStr '.mat']);
save(matFile, 'all_scatterers', 'center_angles', 'aperture_deg', ...
    'fftThreshold', 'fftMinSep', 'fftMaxPeaks', ...
    'lambda_c', 'B', 'c', 'delta_r', 'modelName', '-v7.3');
fprintf('  Data: %s\n', matFile);

% 文本导出
txtFile = fullfile(resultDir, ['fft3d_centers_' nowStr '.txt']);
fid = fopen(txtFile, 'w');
fprintf(fid, '# FFT Peak Detection 3D Scattering Centers\n');
fprintf(fid, '# Model: %s\n', modelName);
fprintf(fid, '# Center angles (deg): %s\n', num2str(center_angles));
fprintf(fid, '# Aperture: %.1f deg, FFT threshold: %.0f dB, Min sep: %.3f m\n', ...
    aperture_deg, fftThreshold, fftMinSep);
fprintf(fid, '# Total centers: %d\n', size(all_scatterers, 1));
fprintf(fid, '#\n# Columns: x(m)  y(m)  z(m)  amplitude  source_theta(deg)\n');
for i = 1:size(all_scatterers, 1)
    fprintf(fid, '  %.6f  %.6f  %.6f  %.6e  %.1f\n', all_scatterers(i, :));
end
fclose(fid);
fprintf('  Centers: %s\n', txtFile);

fprintf('\n========================================\n');
fprintf('  FFT 3D Reconstruction Complete\n');
fprintf('========================================\n');
fprintf('  Model:     %s\n', modelName);
fprintf('  Angles:    %d\n', N_angles);
fprintf('  Centers:   %d\n', size(all_scatterers, 1));
fprintf('  Output:    %s\n', resultDir);
fprintf('========================================\n');

end

%% ========================================================================
%% 辅助函数
%% ========================================================================

function [vertices, faces] = readSTLGeometry(modelFile, ~)
    vertices = []; faces = [];
    if isempty(modelFile) || ~exist(modelFile, 'file')
        fprintf('  STL model not found: %s\n', modelFile);
        return;
    end
    try
        fid = fopen(modelFile, 'rb');
        if fid == -1, return; end
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

function [peaks, I_smoothed] = detectPeaksFFT(I_mag, I_dB, range_axis, ...
    crossrange_axis, threshold_dB, min_sep, max_peaks)

    peaks = zeros(0, 3);
    N_r = length(range_axis);
    N_x = length(crossrange_axis);

    I_smoothed = imgaussfilt(I_dB, 0.6);

    % 3×3 局部最大值
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

    % 亚像素精化
    raw_peaks = zeros(N_raw, 3);
    for k = 1:N_raw
        r_idx = peak_cols(k);
        x_idx = peak_rows(k);

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

    raw_peaks = sortrows(raw_peaks, -3);

    % 合并过近峰
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
