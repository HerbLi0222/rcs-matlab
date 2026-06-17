function total_field = calculate_scattered_field(freq, radius, height, angles_deg)
% 计算雷达回波电场
% 输入：
%   freq        - 雷达工作频率 (Hz)
%   radius      - 圆柱/圆盘半径 (m)
%   height      - 圆柱高度 (m)
%   angles_deg  - 方位角向量 (度)

    % 光速和波长计算
    c = 3e8;
    lambda = c / freq;
    k = 2 * pi / lambda;

    % 角度转为弧度制
    angles_rad = deg2rad(angles_deg);
    angle_count = length(angles_rad);

    % 圆柱反射场计算
    cylinder_field = compute_cylinder_reflection(k, radius, height, angles_rad);

    % 圆盘上下表面反射场计算
    [upper_index, lower_index] = determine_angle_partition(angle_count);
    upper_angles = angles_rad(1:upper_index);
    lower_angles = angles_rad(lower_index:end);

    disk_upper = compute_disk_reflection(k, radius, height, upper_angles);
    disk_lower = compute_disk_reflection(k, radius, height, lower_angles);

    % 合并上下表面反射场
%     size(disk_lower)
%     size(disk_upper)
    disk_total = [disk_upper, disk_lower];

    % 总反射场叠加
    total_field = cylinder_field + disk_total;
end


% 子函数：圆柱体反射场计算
function E = compute_cylinder_reflection(k, a, h, theta)
    term = (1j * k * a * sin(theta) / pi).^0.5;
    phase = exp(2j * k * a * sin(theta));
    E = 0.5 * h * term .* sinc(k * h * cos(theta) / pi) .* phase;
end


% 子函数：圆盘反射场计算
function E = compute_disk_reflection(k, a, h, theta)
    j1 = besselj(1, 2 * k * a * sin(theta));
    phase1=zeros(1,length(theta));
    E=zeros(1,length(theta));
    for i = 1:length(theta)
        if (theta(i)<pi/2)
            phase1(1,i) = exp(-1j * k * h * cos(pi-theta(i)));
        else
            phase1(1,i) = exp(-1j * k * h * cos(theta(i)));
        end
    end
    for i = 1:length(theta)
        if (theta(i)<pi/2)
            E(i) = tan(theta(i)) ;
        else
            E(i) = tan(theta(i)) ;
        end
    end
    E = -1j * 0.5 * a * j1 ./ E .* phase1;
    %E = -1j * 0.5 * a * j1 ./ tan(theta) .* phase1;
    if (theta(1)==0)
        E(1) = -1j * 0.5 * a * (k*a) * exp(-1j * k * h);
    end
end


% 子函数：确定角度分割索引（上半/下半）
function [upper_idx, lower_idx] = determine_angle_partition(N)
    if mod(N, 2) == 0
        upper_idx = N / 2;
    else
        upper_idx = (N + 1) / 2;
    end
    lower_idx = upper_idx + 1;
end
