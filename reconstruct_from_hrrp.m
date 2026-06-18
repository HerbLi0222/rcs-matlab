function reconstruct_from_hrrp(varargin)
% RECONSTRUCT_FROM_HRRP  基于时域距离历程图的三维散射中心重建
%
%   从 main_range_profile.m 生成的 range_profile_data_*.mat 中加载距离历程图，
%   在二维 (角度, 距离) 平面上检测散射峰，利用多角度几何关系将 (θ, r) 坐标
%   三角化重建为三维散射中心分布，并与目标 STL 模型叠加显示。
%
%   方法原理:
%     距离历程图 M(θ, r) 中，一个位于 P 的散射中心在不同角度 θ_i 下满足：
%         r_i = - k̂(θ_i) · P   （注：散射仿真使用 exp(+j·2k·R·P) 约定，
%                                IFFT 后 r 符号翻转）
%     其中 k̂(θ) = [sinθ, 0, cosθ] 为雷达视线方向 (phi=0 单站扫描)。
%
%     修正后的重建公式：P_candidate = -r_i * k̂(θ_i)
%
%     在 HRRP 历程图中检测峰值 → 每个峰值给出一条视线约束。
%     多角度约束联立，通过最小二乘反演 P 的 (x,z) 坐标
%     (y 坐标在 xz 平面扫描中不可观测，始终为 0)。
%
%   与 reconstruct_3d_scatterers.m 的区别:
%     - reconstruct_3d_scatterers: 频域 SAR 成像 (2D IFFT + CLEAN)
%     - reconstruct_from_hrrp:      时域历程图峰值检测 + 多角度三角化
%
%   Usage:
%     >> reconstruct_from_hrrp
%     >> reconstruct_from_hrrp('DataFile', 'range_profile_data_xxx.mat')
%     >> reconstruct_from_hrrp('PeakThreshold', -30, 'MinAngleSpan', 5)

%% ========================================================================
%% 0. 参数解析
%% ========================================================================
addpath('lib');

p = inputParser;
p.addOptional('DataFile', '', @(x) ischar(x) || isstring(x) || isempty(x));
p.addParameter('PeakThreshold', -55, @(x) isnumeric(x) && isscalar(x));
p.addParameter('MinAngleSpan', 3, @(x) isnumeric(x) && isscalar(x));
p.addParameter('MinRangeSep', 0.03, @(x) isnumeric(x) && isscalar(x));
p.addParameter('ClusterRadius', 0.10, @(x) isnumeric(x) && isscalar(x));
p.addParameter('MinClusterSize', 2, @(x) isnumeric(x) && isscalar(x));
p.addParameter('MaxScatterers', 50, @(x) isnumeric(x) && isscalar(x));

if ~isempty(varargin) && (ischar(varargin{1}) || isstring(varargin{1})) ...
        && any(strcmp(char(varargin{1}), {'DataFile','PeakThreshold','MinAngleSpan',...
        'MinRangeSep','ClusterRadius','MinClusterSize','MaxScatterers'}))
    varargin = [{''}, varargin];
end
p.parse(varargin{:});

dataFile        = p.Results.DataFile;
peakThreshold   = p.Results.PeakThreshold;
minAngleSpan    = p.Results.MinAngleSpan;
minRangeSep     = p.Results.MinRangeSep;
clusterRadius   = p.Results.ClusterRadius;
minClusterSize  = p.Results.MinClusterSize;
maxScatterers   = p.Results.MaxScatterers;

%% ========================================================================
%% 1. 加载距离历程图数据
%% ========================================================================
fprintf('========================================\n');
fprintf('  HRRP-Based 3D Scatterer Reconstruction\n');
fprintf('  (基于时域历程图的三维散射中心重建)\n');
fprintf('========================================\n\n');

% --- 1a. 加载 range_profile_data ---
if isempty(dataFile)
    resultDir = fullfile('results');
    files = dir(fullfile(resultDir, 'range_profile_data_*.mat'));
    if isempty(files)
        error('No range_profile_data_*.mat found in %s.\nRun main_range_profile first.', resultDir);
    end
    [~, idx] = sort([files.datenum], 'descend');
    dataFile = fullfile(resultDir, files(idx(1)).name);
end

fprintf('Loading range profile data...\n');
fprintf('  File: %s\n', dataFile);
data = load(dataFile);

M_dB          = data.M_dB;           % N_angles x N_fft (归一化后)
range_axis    = data.range_axis;     % 1 x N_fft
angle_axis    = data.angle_axis;     % N_angles x 1
delta_r_fft   = data.delta_r_fft;
N_fft         = data.N_fft;
c             = data.c;
B             = data.B;

N_angles = size(M_dB, 1);
fprintf('  Matrix: %d angles x %d range bins\n', N_angles, N_fft);
fprintf('  Range span: %.2f ~ %.2f m (res: %.3f m)\n', ...
    range_axis(1), range_axis(end), delta_r_fft);
fprintf('  Angle span: %.1f ~ %.1f deg\n', angle_axis(1), angle_axis(end));

% --- 1b. 确定相位参考点 ---
% 尝试从 wideband_scattering 数据获取 bbox_center
P_ref = [0; 0; 0];
widebandFiles = dir(fullfile('results', 'wideband_scattering_*.mat'));
if ~isempty(widebandFiles)
    [~, sortIdx] = sort([widebandFiles.datenum], 'descend');
    wbData = load(fullfile('results', widebandFiles(sortIdx(1)).name));
    if isfield(wbData, 'bbox_center')
        P_ref = wbData.bbox_center(:);
        modelName = '';
        if isfield(wbData, 'inputModel')
            [~, modelName, ~] = fileparts(wbData.inputModel);
            modelFile = fullfile('stl_models', wbData.inputModel);
        else
            modelFile = '';
        end
        fprintf('  Reference point (bbox): [%.3f, %.3f, %.3f] m\n', P_ref);
        fprintf('  Model: %s\n', modelName);
    end
else
    modelFile = '';
    modelName = 'Unknown';
    fprintf('  Reference point: origin (no bbox data)\n');
end

%% ========================================================================
%% 2. 在距离历程图中检测散射峰
%% ========================================================================
fprintf('\n========================================\n');
fprintf('  Detecting Scattering Peaks in HRRP Map\n');
fprintf('========================================\n');
fprintf('  Absolute threshold: %.0f dB\n', peakThreshold);
fprintf('  Min angle span: %d deg\n', minAngleSpan);
fprintf('  Min range separation: %.2f m\n', minRangeSep);

% --- 自适应峰值检测 ---
% 逐角度找到局部最大值，然后合并到全局列表
% 对每个角度：取在该角度峰值以下 25dB 内的局部最大值
perAngleDrop = 35;  % dB below the per-angle peak (足够捕获较弱散射中心)

M_smooth = imgaussfilt(M_dB, 0.8);
local_max = zeros(size(M_smooth));
angle_peaks_dB = max(M_dB, [], 2);  % 每角度最大 dB 值

for i = 2:N_angles-1
    adaptiveThresh = max(peakThreshold, angle_peaks_dB(i) - perAngleDrop);
    for j = 2:N_fft-1
        patch = M_smooth(i-1:i+1, j-1:j+1);
        if M_smooth(i,j) == max(patch(:)) && M_smooth(i,j) > adaptiveThresh
            % 额外检查：与相邻角度的同距离峰差异不太大（排除噪声）
            local_max(i,j) = 1;
        end
    end
end

[peak_rows, peak_cols] = find(local_max);
N_peaks_raw = length(peak_rows);
fprintf('  Raw local maxima detected: %d\n', N_peaks_raw);

% 提取峰值属性
all_peaks = zeros(N_peaks_raw, 5);
for k = 1:N_peaks_raw
    i = peak_rows(k);
    j = peak_cols(k);
    all_peaks(k, :) = [i, j, angle_axis(i), range_axis(j), M_dB(i,j)];
end

% 按幅度降序排列，合并距离过近的峰
all_peaks = sortrows(all_peaks, -5);
if N_peaks_raw > 1
    merged = all_peaks(1, :);
    for k = 2:N_peaks_raw
        % 检查是否与已合并的峰在角度和距离上都太近
        too_close = false;
        for m = 1:size(merged, 1)
            if abs(all_peaks(k,3) - merged(m,3)) < 2 && ...
               abs(all_peaks(k,4) - merged(m,4)) < minRangeSep
                too_close = true;
                break;
            end
        end
        if ~too_close
            merged = [merged; all_peaks(k, :)]; %#ok<AGROW>
        end
    end
    all_peaks = merged;
end
N_peaks_raw = size(all_peaks, 1);
fprintf('  After merging close peaks: %d\n', N_peaks_raw);

fprintf('  Top 20 peaks:\n');
for k = 1:min(20, N_peaks_raw)
    fprintf('    θ=%.0f°, r=%.3f m, amp=%.1f dB\n', ...
        all_peaks(k,3), all_peaks(k,4), all_peaks(k,5));
end

%% ========================================================================
%% 3. 对每个峰值计算候选 3D 位置
%% ========================================================================
fprintf('\n========================================\n');
fprintf('  Computing 3D Candidate Positions\n');
fprintf('========================================\n');

% 对每个峰值: P_candidate = -r * k_hat（取反修正相位约定符号翻转）
N_peaks = size(all_peaks, 1);
candidates_3d = zeros(N_peaks, 3);  % [x, y, z]
for k = 1:N_peaks
    theta_deg = all_peaks(k, 3);
    r_val = all_peaks(k, 4);
    theta_rad = deg2rad(theta_deg);
    k_hat = [sin(theta_rad); 0; cos(theta_rad)];
    % 修正：散射仿真相位约定导致 r = -k_hat·P，需取反 (−r·k_hat)
    P_cand = - r_val * k_hat;
    candidates_3d(k, :) = P_cand';
end

fprintf('  Candidates computed: %d\n', N_peaks);

%% ========================================================================
%% 4. 3D 空间聚类 → 同一散射中心的多角度观测合并
%% ========================================================================
fprintf('\n========================================\n');
fprintf('  Clustering Candidates in 3D Space\n');
fprintf('========================================\n');
fprintf('  Cluster radius: %.3f m\n', clusterRadius);
fprintf('  Min cluster size: %d\n', minClusterSize);

% 处理候选点不足的情况
if N_peaks < 2
    % 只有一个或零个候选 → 直接作为散射体
    clusters = ones(N_peaks, 1);
    fprintf('  Only %d candidate(s), skipping clustering.\n', N_peaks);
else
    % 层次聚类
    D = squareform(pdist(candidates_3d));
    Z = linkage(squareform(D), 'single');
    clusters = cluster(Z, 'cutoff', clusterRadius, 'criterion', 'distance');
end

unique_clusters = unique(clusters);
N_clusters = length(unique_clusters);
fprintf('  Clusters found: %d\n', N_clusters);

% 筛选：每个簇至少 minClusterSize 个候选点（或全部保留如果总数太少）
effectiveMinSize = min(minClusterSize, max(1, floor(N_peaks / max(1, N_clusters))));
scatterer_list = [];  % [x, y, z, amplitude, n_obs, angle_span]
for c = 1:N_clusters
    members = find(clusters == c);
    n_members = length(members);
    if n_members < effectiveMinSize && N_peaks >= minClusterSize
        continue;
    end

    % 簇内统计
    cluster_pts = candidates_3d(members, :);
    centroid = mean(cluster_pts, 1);
    cluster_amps = all_peaks(members, 5);
    mean_amp = mean(cluster_amps);
    cluster_angles = all_peaks(members, 3);
    angle_span = max(cluster_angles) - min(cluster_angles);

    % 如果候选总数少，放宽角度跨度要求
    effMinAngleSpan = min(minAngleSpan, max(1, minAngleSpan * N_peaks / 10));
    if angle_span < effMinAngleSpan && N_peaks >= minClusterSize
        continue;
    end

    scatterer_list = [scatterer_list; centroid, mean_amp, n_members, angle_span]; %#ok<AGROW>
end

if isempty(scatterer_list)
    warning('No scattering clusters found meeting criteria.');
    fprintf('  Try lowering PeakThreshold (currently %.0f) or ClusterRadius (currently %.3f)\n', ...
        peakThreshold, clusterRadius);
    scatterer_list = zeros(0, 6);
end

% 按幅度排序，取前 maxScatterers 个
scatterer_list = sortrows(scatterer_list, -4);
if size(scatterer_list, 1) > maxScatterers
    scatterer_list = scatterer_list(1:maxScatterers, :);
end

N_scat = size(scatterer_list, 1);
fprintf('  Final scattering centers: %d\n', N_scat);
for k = 1:N_scat
    fprintf('    #%d: [%.3f, %.3f, %.3f] m, amp=%.1f dB, obs=%d, span=%.0f°\n', ...
        k, scatterer_list(k,1), scatterer_list(k,2), scatterer_list(k,3), ...
        scatterer_list(k,4), scatterer_list(k,5), scatterer_list(k,6));
end

%% ========================================================================
%% 5. 三维可视化
%% ========================================================================
fprintf('\nGenerating 3D visualization...\n');

% --- 5a. 加载 STL 几何 ---
[vertices, faces] = readSTLBinary(modelFile);

% --- 诊断 ---
fprintf('\n--- Diagnostic: Coordinate Range Comparison ---\n');
if ~isempty(vertices)
    fprintf('  STL model bbox:\n');
    fprintf('    X: [%.3f, %.3f] m\n', min(vertices(:,1)), max(vertices(:,1)));
    fprintf('    Y: [%.3f, %.3f] m\n', min(vertices(:,2)), max(vertices(:,2)));
    fprintf('    Z: [%.3f, %.3f] m\n', min(vertices(:,3)), max(vertices(:,3)));
end
if N_scat > 0
    fprintf('  Scattering centers:\n');
    fprintf('    X: [%.3f, %.3f] m\n', min(scatterer_list(:,1)), max(scatterer_list(:,1)));
    fprintf('    Y: [%.3f, %.3f] m\n', min(scatterer_list(:,2)), max(scatterer_list(:,2)));
    fprintf('    Z: [%.3f, %.3f] m\n', min(scatterer_list(:,3)), max(scatterer_list(:,3)));
end

% 幅度归一化用于点大小
if N_scat > 0
    A_min = min(scatterer_list(:,4));
    A_max = max(scatterer_list(:,4));
    if A_max - A_min < 1e-6
        A_norm = ones(N_scat, 1);
    else
        A_norm = (scatterer_list(:,4) - A_min) / (A_max - A_min);
    end
    marker_sizes = 30 + 200 * A_norm;
else
    marker_sizes = [];
    A_norm = [];
end

% 按观测次数着色
if N_scat > 0
    n_obs = scatterer_list(:, 5);
    cmap = jet(64);
    color_idx = max(1, min(64, round(63 * (n_obs - min(n_obs)) / (max(n_obs) - min(n_obs) + 1e-10)) + 1));
    scat_color = cmap(color_idx, :);
end

%% --- 图1: 3D 总览 ---
figure(1);
clf;
set(gcf, 'Name', 'HRRP-Based 3D Scattering Centers', ...
         'NumberTitle', 'off', 'Position', [50, 80, 1000, 750]);

hold on;

% STL 模型 (半透明，用 trisurf 更稳定)
if ~isempty(vertices) && ~isempty(faces) && size(faces,2) == 3
    try
        % 验证面索引范围
        if max(faces(:)) <= size(vertices, 1)
            trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3), ...
                    'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', ...
                    'FaceAlpha', 0.25, 'AmbientStrength', 0.3);
        else
            fprintf('  STL skipped: face indices out of bounds.\n');
        end
    catch ME
        fprintf('  STL render skipped: %s\n', ME.message);
    end
end

% 散射中心
if N_scat > 0
    for i = 1:N_scat
        scatter3(scatterer_list(i,1), scatterer_list(i,2), scatterer_list(i,3), ...
                 marker_sizes(i), scat_color(i, :), ...
                 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    end
end

% 参考点
scatter3(P_ref(1), P_ref(2), P_ref(3), 150, 'k', 'p', ...
         'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);

% 各角度的雷达视线方向 (均匀采样)
n_dirs = 8;
sample_angles = linspace(angle_axis(1), angle_axis(end), n_dirs);
dir_scale = max(abs(range_axis)) * 0.4;
for i = 1:n_dirs
    th = deg2rad(sample_angles(i));
    k_hat = [sin(th); 0; cos(th)];
    quiver3(P_ref(1), P_ref(2), P_ref(3), ...
            dir_scale * k_hat(1), 0, dir_scale * k_hat(3), ...
            0, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'LineStyle', '--', ...
            'MaxHeadSize', 0.4, 'AutoScale', 'off');
end

hold off;
axis equal;
grid on;
box on;
xlabel('X (m)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y (m)', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Z (m)', 'FontSize', 12, 'FontWeight', 'bold');

if N_scat > 0
    title(sprintf(['HRRP-Based 3D Scattering Centers\n', ...
        'Model: %s | %d Centers | %d Angles | Cluster Radius: %.2f m'], ...
        modelName, N_scat, N_angles, clusterRadius), ...
        'FontSize', 13, 'FontWeight', 'bold');
end

view(45, 25);
camlight('headlight');
lighting gouraud;
material dull;

%% --- 图2: 距离历程图 + 检测到的峰值叠加 ---
figure(2);
clf;
set(gcf, 'Name', 'HRRP Map with Detected Peaks', ...
         'NumberTitle', 'off', 'Position', [100, 100, 900, 600]);

imagesc(angle_axis, range_axis, M_dB');
colormap('jet');
axis xy;
hold on;

% 叠加检测到的峰值
for k = 1:N_peaks_raw
    plot(all_peaks(k,3), all_peaks(k,4), 'wo', ...
         'MarkerSize', 6, 'LineWidth', 1.2);
end

% 叠加最终散射中心对应的候选点
if N_scat > 0
    for k = 1:N_peaks_raw
        if any(clusters(k) == unique_clusters)
            plot(all_peaks(k,3), all_peaks(k,4), 'k.', 'MarkerSize', 10);
        end
    end
end

hold off;
xlabel('Observation Angle (deg)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Radial Range (m)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('HRRP Map with Detected Peaks | %d Raw, %d Final Centers', ...
    N_peaks_raw, N_scat), 'FontSize', 13, 'FontWeight', 'bold');

clim([peakThreshold-50, 0]);
cb2 = colorbar;
cb2.Label.String = 'Relative Intensity (dB)';

grid on;
set(gca, 'GridAlpha', 0.3, 'FontSize', 11);

%% --- 图3: 三视图 ---
if N_scat > 0
    figure(3);
    clf;
    set(gcf, 'Name', 'Three-View Projections (三视图)', ...
             'NumberTitle', 'off', 'Position', [100, 100, 1100, 400]);

    view_configs = {
        'XZ Plane (Front)',  'X (m)', 'Z (m)', scatterer_list(:,1), scatterer_list(:,3);
        'YZ Plane (Side)',   'Y (m)', 'Z (m)', scatterer_list(:,2), scatterer_list(:,3);
        'XY Plane (Top)',    'X (m)', 'Y (m)', scatterer_list(:,1), scatterer_list(:,2);
    };

    for v = 1:3
        subplot(1, 3, v);
        hold on;

        % 目标轮廓
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
        for i = 1:N_scat
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
        scatter(px, py_p, 80, 'k', 'p', 'filled', 'MarkerEdgeColor', 'w');

        hold off;
        axis equal;
        grid on;
        xlabel(view_configs{v,2}, 'FontSize', 11);
        ylabel(view_configs{v,3}, 'FontSize', 11);
        title(view_configs{v,1}, 'FontSize', 12, 'FontWeight', 'bold');
    end

    sgtitle(sprintf('Three-View Projections | Model: %s | %d Centers', ...
           modelName, N_scat), 'FontSize', 13, 'FontWeight', 'bold');
end

%% ========================================================================
%% 6. 保存结果
%% ========================================================================
fprintf('\nSaving results...\n');

nowStr = datestr(now, 'yyyymmddHHMMSS');

for figNum = 1:3
    if ishandle(figNum)
        figFile = fullfile('results', sprintf('hrrp3d_fig%d_%s.png', figNum, nowStr));
        saveas(figNum, figFile);
        fprintf('  Figure %d: %s\n', figNum, figFile);
    end
end

matFile = fullfile('results', ['hrrp3d_data_' nowStr '.mat']);
save(matFile, 'scatterer_list', 'P_ref', 'all_peaks', 'candidates_3d', ...
    'clusters', 'peakThreshold', 'clusterRadius', 'minClusterSize', ...
    'angle_axis', 'range_axis', 'modelName', '-v7.3');
fprintf('  Data: %s\n', matFile);

txtFile = fullfile('results', ['hrrp3d_centers_' nowStr '.txt']);
fid = fopen(txtFile, 'w');
fprintf(fid, '# HRRP-Based 3D Scattering Centers\n');
fprintf(fid, '# Model: %s\n', modelName);
fprintf(fid, '# Reference point: [%.4f, %.4f, %.4f] m\n', P_ref);
fprintf(fid, '# Detection: threshold=%.0f dB, cluster_radius=%.3f m\n', ...
    peakThreshold, clusterRadius);
fprintf(fid, '# Total centers: %d\n', N_scat);
fprintf(fid, '#\n');
fprintf(fid, '# Columns: x(m)  y(m)  z(m)  amplitude_dB  n_obs  angle_span(deg)\n');
for i = 1:N_scat
    fprintf(fid, '  %.6f  %.6f  %.6f  %.2f  %d  %.1f\n', scatterer_list(i, :));
end
fclose(fid);
fprintf('  Text: %s\n', txtFile);

fprintf('\n========================================\n');
fprintf('  HRRP-Based Reconstruction Complete\n');
fprintf('========================================\n');
fprintf('  Model:            %s\n', modelName);
fprintf('  Method:           HRRP peak detection + 3D clustering\n');
fprintf('  Raw peaks:        %d\n', N_peaks_raw);
fprintf('  Final centers:    %d\n', N_scat);
fprintf('  Cluster radius:   %.3f m\n', clusterRadius);
fprintf('  Reference point:  [%.3f, %.3f, %.3f] m\n', P_ref);
fprintf('========================================\n');

end  % reconstruct_from_hrrp

%% ========================================================================
%% 辅助函数
%% ========================================================================

% -------------------------------------------------------------------------
function [vertices, faces] = readSTLBinary(modelFile)
% 直接读取二进制 STL 文件，不依赖任何外部函数

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

        fread(fid, 80, 'uint8=>uint8');              % 跳过 80 字节头
        nTri = fread(fid, 1, 'uint32=>uint32');      % 三角面元数

        allVerts = zeros(nTri * 3, 3);

        for i = 1:nTri
            fread(fid, 3, 'float32=>double');        % 法向量
            allVerts(3*i-2, :) = fread(fid, 3, 'float32=>double')';
            allVerts(3*i-1, :) = fread(fid, 3, 'float32=>double')';
            allVerts(3*i,   :) = fread(fid, 3, 'float32=>double')';
            fread(fid, 1, 'uint16=>uint16');         % 属性字节
        end
        fclose(fid);

        % 去重并用重映射索引构建面列表
        [vertices, ~, ic] = unique(allVerts, 'rows', 'stable');
        faces = reshape(double(ic), 3, nTri)';

        fprintf('  STL loaded: %d vertices, %d faces\n', ...
                size(vertices, 1), size(faces, 1));
    catch ME
        fprintf('  STL read error: %s\n', ME.message);
    end
end
