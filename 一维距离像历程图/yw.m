clear,clc;

% 光速定义
c = 3e8;

% 雷达系统参数
center_freq = 3e9;
bandwidth = 1e9;
low_freq = center_freq - bandwidth / 2;
high_freq = center_freq + bandwidth / 2;
freq_points = 100;
angle_resolution = 0.01;
azimuth_angles = 0:angle_resolution:180;
angle_count = length(azimuth_angles);

% 几何体配置
cylinder_radius = 2;
cylinder_height = 4*3^0.5 ;

% 距离像相关设置
oversampling_factor = 4;
distance_step = c / (2 * bandwidth) / oversampling_factor;
range_profile = (1:oversampling_factor * freq_points) * distance_step - (oversampling_factor * freq_points / 2) * distance_step;

% 频率点生成
frequencies = linspace(low_freq, high_freq, freq_points);

% 回波信号初始化
echo_signal = zeros(freq_points, angle_count);
% 一维距离像初始化
echo_fft = zeros(angle_count, oversampling_factor * freq_points);

for idx_freq = 1:freq_points
    % 计算回波电场
    echo_signal(idx_freq, :) = calculate_scattered_field(frequencies(idx_freq), cylinder_radius, cylinder_height, azimuth_angles);
end

for idx_theta = 1:angle_count
    % 不同角度下的距离像计算
    echo_fft(idx_theta,:) = ifftshift(ifft(echo_signal(:,idx_theta)',oversampling_factor * freq_points));
end

% 指定角度处理
radar_angle = 30;
angle_index = round(radar_angle * (angle_count - 1) / 180) + 1;

% 绘制回波电场幅度频谱
figure(1);
plot(frequencies/1e9,abs(echo_signal(:,angle_index)'));
xlabel('f (GHz)');                   % x轴标签
ylabel('Es (V/m)');                  % y轴标签
title('回波电场幅度频谱');            % 图标题 

selected_fft = echo_fft(angle_index,:)';
range_profile_data = abs(selected_fft);
% 绘图：一维距离像
figure(2);
range_profile_data_db = 20 * log10(range_profile_data);
plot(range_profile, range_profile_data_db); 
xlim([min(range_profile), max(range_profile)]);
xlabel('Distance (m)');
ylabel('Amplitude (dB)');
title('1D Range Profile,\theta=0°');


echo_fft_2d = abs(echo_fft)';

% 创建图像数据范围
x = azimuth_angles;          % 横轴：角度
y = range_profile;           % 纵轴：距离

% 使用 imagesc 或 pcolor 绘图
figure(3);
imagesc(x, y, 20*log10(echo_fft_2d)); % 显示图像
axis xy;                     % Y轴正方向朝上
colormap('jet');             % 设置颜色映射（红亮蓝暗）
colorbar;                    % 显示颜色条

% 添加标签与标题
xlabel('Azimuth Angle (°)');
ylabel('Range (m)');
title('2D Range-Angle Profile (Echo Intensity in dB)');

% 可选：设置坐标轴范围
xlim([min(x), max(x)]);
ylim([min(y), max(y)]);

% 可选：添加网格线
grid on;

% 如果你想让蓝色代表小值、红色代表大值，也可以自定义 colormap：
% cmap = flipud(colormap('jet'));  % 翻转颜色映射
% colormap(cmap);

