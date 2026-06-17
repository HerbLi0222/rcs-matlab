clear, clc;
% 光速定义
c = 3e8;
% 雷达系统参数
center_freq = 3e9;
bandwidth =1e9;
freq_points = 101;
fre_interval = bandwidth/(freq_points-1);
angle_resolution = 0.01;
start_theta = 0; % 角度为度数
end_theta = 180 ;   % 角度为度数
% 频率点生成
low_freq = center_freq - bandwidth / 2;
high_freq = center_freq + bandwidth / 2;
frequencies = linspace(low_freq, high_freq, freq_points);
%角度生成
azimuth_angles = start_theta:angle_resolution:end_theta; % 度数
angle_count = length(azimuth_angles);
dx = c/bandwidth/2;
x_width = c/2/fre_interval;
dy = c/center_freq/2/deg2rad(end_theta-start_theta);
y_width = c/2/center_freq/deg2rad(angle_resolution);
% 几何体配置
cylinder_radius = 2;
cylinder_height = 4*3^0.5;

% --- 生成原始 kx 和 ky 坐标以及回波信号 ---
[AZ_ANGLES_GRID, FREQ_GRID] = meshgrid(deg2rad(azimuth_angles), frequencies);

Kx_original = 4*pi/c * FREQ_GRID .* cos(AZ_ANGLES_GRID); 
Ky_original = 4*pi/c * FREQ_GRID .* sin(AZ_ANGLES_GRID);

echo_signal = zeros(freq_points, angle_count);

for idx_freq = 1:freq_points
    % 替换为你的实际 calculate_scattered_field 函数调用
    echo_signal(idx_freq, :) = calculate_scattered_field(frequencies(idx_freq), cylinder_radius, cylinder_height, azimuth_angles);
end

% 将 Kx_original, Ky_original, echo_signal 展平为列向量，作为散点数据
points_kx = Kx_original(:);
points_ky = Ky_original(:);
values_echo = echo_signal(:);

% --- 插值变换：使用 scatteredInterpolant ---
% 1. 创建 scatteredInterpolant 对象
% 'linear' 表示使用线性插值，也可以尝试 'nearest', 'natural', 'cubic'
F = scatteredInterpolant(points_kx, points_ky, values_echo, 'linear', 'none');

% 2. 定义新的均匀网格的范围和点数
kx_min = min(points_kx);
kx_max = max(points_kx);
ky_min = min(points_ky);
ky_max = max(points_ky);

% 决定插值后网格的点数。可以根据原始数据尺寸或期望分辨率设定。
% 为了保持与原始数据差不多的分辨率，可以继续使用之前的方法：
N_interp = angle_count; % 新网格的列数 (kx方向)
M_interp = freq_points;  % 新网格的行数 (ky方向)

% 3. 生成新的均匀 kx 和 ky 轴
k_x_new_axis = linspace(kx_min, kx_max, N_interp);
k_y_new_axis = linspace(ky_min, ky_max, M_interp);

% 4. 使用 meshgrid 生成新的均匀查询网格
[K_X_NEW_QUERY, K_Y_NEW_QUERY] = meshgrid(k_x_new_axis, k_y_new_axis);

% 5. 使用 scatteredInterpolant 对象在新的查询网格上进行插值
echo_signal_interp = F(K_X_NEW_QUERY, K_Y_NEW_QUERY);

% 7. 根据图片中“新网格未插值的地方是取0”的规则，将 NaN 值替换为 0
echo_signal_interp(isnan(echo_signal_interp)) = 0;

%绘制散射场热力图
figure;
imagesc(20*log10(abs(echo_signal_interp')));  % 使用绝对值显示强度
colorbar;                    % 添加颜色条
xlabel('频率索引', 'fontsize', 12);
ylabel('角度索引', 'fontsize', 12);
title('散射场强度热力图（频率 vs 角度）', 'fontsize', 14);
colormap('jet');             % 可选：使用 jet 颜色映射
axis xy;                     % 确保第1行对应最低频率

% 二维傅里叶变换
isar_image = fftshift(fft2(echo_signal_interp')); %angle x fre

% 取绝对值并转换为dB单位
isar_image = abs(isar_image);
isar_image_db = 20 * log10(isar_image); 

% 显示ISAR图像
figure;
ff=linspace(-(freq_points-1)/2*dx,(freq_points-1)/2*dx,freq_points);
pp=linspace(-(angle_count-1)/2*dy,(angle_count-1)/2*dy,angle_count);
imagesc(ff, pp, isar_image_db);
axis xy;
colormap('jet');
colorbar;
axis([-x_width/2 x_width/2 -y_width/2 y_width/2]);
axis square
xlabel('距离向 (m)');
ylabel('方位向 (m)');
title('ISAR Image (Intensity in dB)');
% clim = get(gca,'CLim');
% set(gca,'CLim',clim(2) + [-20 0]);