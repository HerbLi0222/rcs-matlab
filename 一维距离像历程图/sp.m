clear,clc;

% 光速定义
c = 3e8;

% 雷达系统参数
center_freq = 3e9;
angle_resolution = 1;
N=2^11;
azimuth_angles = linspace(0,180,N);

% 几何体配置
cylinder_radius = 2;

cylinder_height = 4*3^0.5;

% 回波信号初始化
echo_signal = zeros(1, N);
echo_signal(1, :) = calculate_scattered_field(center_freq, cylinder_radius, cylinder_height, azimuth_angles);

t = 1:N;
[tfr, ~, ~] = tfrspwv(echo_signal',t,N);
reference = 20*log10(abs(ifftshift(tfr,1)));
w=100;
fD_min = -w/(2*angle_resolution);
fD_max = w/(2*angle_resolution);
imagesc([0,180],[fD_min,fD_max],reference);
colorbar();
colormap(jet);
title(sprintf('a=2,h=4√3 半球体-圆柱时频像\n 入射频率为 3GHz'));
xlabel('θ /°');
ylabel('Doppler frequency / Hz');
