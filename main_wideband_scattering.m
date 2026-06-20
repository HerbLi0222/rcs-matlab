% MAIN_WIDEBAND_SCATTERING  宽带复散射场计算程序
%
%   Wideband Complex Scattering Field Calculation using Physical Optics (PO)
%
%   根据"距离历程计算方法.txt"第2节的流程，对每个观测角度和每个频率点
%   计算远区复数散射场 E_s(f_j, θ_i)，保留幅度和相位信息。
%
%   计算流程:
%     1. 读取宽带仿真参数（频率范围、频点数、角度范围、角度步长等）
%     2. 提取目标几何数据（STL三角面元模型）
%     3. 加载材料属性（可选）
%     4. 预计算几何量（与频率无关）
%     5. 双循环：扫角 → 扫频 → 计算复散射场（含相位、遮挡判断）
%     6. 保存复数散射矩阵 S_complex(θ, f) 到 results/ 目录
%
%   输出文件:
%     results/wideband_scattering_<timestamp>.mat
%       包含: S_complex (复散射场矩阵), freq_array, theta_array, phi_array,
%             param (参数摘要), wave_array, input_model
%
%   Usage:
%     >> main_wideband_scattering

clear; clc;

% Add library path
addpath('lib');

fprintf('========================================\n');
fprintf('  Wideband Scattering Field Calculation\n');
fprintf('  (宽带复散射场计算)\n');
fprintf('========================================\n\n');

%% 1. 加载宽带仿真参数
fprintf('Loading wideband parameters...\n');

inputDataFile = fullfile('input_files', 'input_data_file_wideband.txt');
fid = fopen(inputDataFile, 'r');
if fid == -1
    error('main_wideband_scattering:FileNotFound', ...
          'Cannot open: %s', inputDataFile);
end

% 按行读取参数（与 getParamsFromFile 相同的解析逻辑）
paramRaw = {};
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if isempty(line) || line(1) == '#'
        continue;
    end
    [val, status] = str2num(line); %#ok<ST2NM>
    if status && ~isempty(val)
        paramRaw{end+1} = val; %#ok<AGROW>
    else
        paramRaw{end+1} = line; %#ok<AGROW>
    end
end
fclose(fid);

% 解析参数
inputModel  = paramRaw{1};       % STL 模型文件名
freq_start  = paramRaw{2} * 1e9; % 起始频率 (GHz -> Hz)
freq_stop   = paramRaw{3} * 1e9; % 终止频率 (GHz -> Hz)
N_f         = paramRaw{4};       % 频率点数
corr        = paramRaw{5};       % 相关距离 (m)
delstd      = paramRaw{6};       % 表面粗糙度标准差 (m)
ipol        = paramRaw{7};       % 极化方式 (0=TM-z, 1=TE-z)
rs          = paramRaw{8};       % 材料类型 (0=PEC, 1=material specific)
pstart      = paramRaw{9};       % 起始 phi 角 (deg)
pstop       = paramRaw{10};      % 终止 phi 角 (deg)
delp        = paramRaw{11};      % phi 步长 (deg)
tstart      = paramRaw{12};      % 起始 theta 角 (deg)
tstop       = paramRaw{13};      % 终止 theta 角 (deg)
delt        = paramRaw{14};      % theta 步长 (deg)
matrlpath   = paramRaw{end};     % 材料文件路径

MATERIALESPECIFICO = 1;

% 验证频点数（建议为 2 的幂，便于 FFT）
if N_f > 0 && (bitand(N_f, N_f - 1) ~= 0)
    fprintf('  Warning: N_f=%d is not a power of 2. Consider using 128, 256, 512, etc.\n', N_f);
end

% 生成频率数组
freq_array = linspace(freq_start, freq_stop, N_f);
freq_step = (freq_stop - freq_start) / (N_f - 1);
B = freq_stop - freq_start;  % 带宽

fprintf('Parameters loaded for wideband simulation.\n');
fprintf('  Model:       %s\n', inputModel);
fprintf('  Frequency:   %.2f ～ %.2f GHz (%d points)\n', freq_start/1e9, freq_stop/1e9, N_f);
fprintf('  Bandwidth:   %.2f GHz\n', B/1e9);
fprintf('  Polarization: %s\n', char('TM-z' + (ipol==1)*('TE-z'-'TM-z')));
fprintf('  Phi:    %.1f:%.1f:%.1f deg\n', pstart, delp, pstop);
fprintf('  Theta:  %.1f:%.1f:%.1f deg\n', tstart, delt, tstop);

%% 2. 提取目标几何数据
fprintf('Processing geometry...\n');

% 转换 STL 模型为坐标和面元数据
modelFile = fullfile('stl_models', inputModel);
stlConverter(modelFile);

coordinatesData = extractCoordinatesData(rs);

x     = coordinatesData{1};
y     = coordinatesData{2};
z     = coordinatesData{3};
xpts  = coordinatesData{4};
ypts  = coordinatesData{5};
zpts  = coordinatesData{6};
nverts = coordinatesData{7};
nfc   = coordinatesData{8};
node1 = coordinatesData{9};
node2 = coordinatesData{10};
node3 = coordinatesData{11};
iflag = coordinatesData{12};
ilum  = coordinatesData{13};
Rs    = coordinatesData{14};
ntria = coordinatesData{15};
vind  = coordinatesData{16};
r     = coordinatesData{17};

% ---- 计算包围盒中心（对任意多面体通用，不受网格密度分布影响）----
% 算术平均（顶点质心）在网格不均匀时会有偏——包围盒中心更鲁棒
bbox_center = [(min(x)+max(x))/2, (min(y)+max(y))/2, (min(z)+max(z))/2];
fprintf('  Bounding box center: (%.3f, %.3f, %.3f) m\n', bbox_center);
fprintf('  Vertex centroid:     (%.3f, %.3f, %.3f) m\n', mean(x), mean(y), mean(z));

fprintf('  Triangles: %d\n', ntria);

%% 3. 加载材料属性
matrl = {};
if rs == MATERIALESPECIFICO
    try
        matrl = getEntrysFromMatrlFile(ntria, matrlpath);
        fprintf('Material file loaded: %s\n', matrlpath);
    catch ME
        error('Material file error: %s', ME.message);
    end
end

%% 4. 物理常数与预计算
c = 3e8;  % 光速 (m/s)
rad = pi / 180;  % 度转弧度

% 距离分辨率检查（必须满足 L > 15 * Δr，即 B > 15*c/(2*L)）
delta_r = c / (2 * B);
fprintf('\nRange resolution check:\n');
fprintf('  Bandwidth B = %.2f GHz\n', B/1e9);
fprintf('  Range resolution Δr = c/(2B) = %.4f m\n', delta_r);
fprintf('  Required: L > 15*Δr = %.4f m (ensure target fits in range window)\n', 15*delta_r);

% 预计算几何量（与频率无关）
fprintf('Pre-computing geometry arrays...\n');
[Area, alpha, beta, N, d, ip, it] = calculateValues(pstart, pstop, delp, ...
    tstart, tstop, delt, ntria, rad);
[N, d, Area, beta, alpha] = productVector(ntria, N, r, d, Area, alpha, beta, vind);

% 角度计数器
if ip == 0, ip = 1; end
if it == 0, it = 1; end
N_angles_total = ip * it;  % 总观测角度数

fprintf('  Angular grid: %d phi × %d theta = %d obs. angles\n', ip, it, N_angles_total);

% 预分配：复数散射场矩阵 S_complex (N_angles × N_f)
% 每行对应一个观测角度，每列对应一个频点
S_complex = zeros(N_angles_total, N_f);

% 存储角度数组
phi_array   = zeros(ip, it);
theta_array = zeros(ip, it);

%% 5. 预计算入射场极化（与角度有关但频率无关的部分）
[pol, Et, Ep] = getPolarization(ipol);
Co = 1;  % 顶点波幅

%% 6. 宽带复散射场计算 —— 双循环：扫角 → 扫频
fprintf('\n========================================\n');
fprintf('  Computing Wideband Scattering Fields\n');
fprintf('========================================\n');
fprintf('  Total: %d angles × %d frequencies = %d evaluations\n', ...
    N_angles_total, N_f, N_angles_total * N_f);
tic;

% 角度编号（扁平化索引）
i_angle = 0;

% 进度显示控制
progress_interval = max(1, floor(N_angles_total / 20));
next_progress     = progress_interval;
start_time        = tic;

for i1 = 1:ip      % ---- phi 角循环 ----
    for i2 = 1:it  % ---- theta 角循环 ----

        i_angle = i_angle + 1;

        % ----- 观测角度 -----
        if ip > 1
            phi_val   = pstart + (i1 - 1) * delp;
        else
            phi_val   = pstart;
        end
        if it > 1
            theta_val = tstart + (i2 - 1) * delt;
        else
            theta_val = tstart;
        end

        phi_array(i1, i2)   = phi_val;
        theta_array(i1, i2) = theta_val;

        phr = phi_val * rad;
        thr = theta_val * rad;

        % ----- 进度显示（约 5% 步进）-----
        if i_angle >= next_progress || i_angle == 1 || i_angle == N_angles_total
            elapsed  = toc(start_time);
            pct      = 100 * i_angle / N_angles_total;
            rate     = i_angle / elapsed;
            eta      = (N_angles_total - i_angle) / max(rate, 1e-6);
            fprintf('  [%3.0f%%] angle %d/%d  theta=%.1f  phi=%.1f  |  elapsed: %.1fs  ETA: %.1fs\n', ...
                pct, i_angle, N_angles_total, theta_val, phi_val, elapsed, eta);
            next_progress = next_progress + progress_interval;
        end

        % ----- 全局方向余弦（与频率无关）-----
        % globalAngles 需要数组输入进行索引赋值，创建局部哑数组
        U_dummy = zeros(1, 1);
        V_dummy = zeros(1, 1);
        W_dummy = zeros(1, 1);
        [~, ~, ~, D0, uu, vv, ww, u, v, w] = ...
            globalAngles(U_dummy, V_dummy, W_dummy, thr, phr, 1, 1);

        % 雷达视线方向单位矢量
        R_vec = [u; v; w];

        % ----- 入射场在全局笛卡尔坐标（与频率无关）-----
        % incidentFieldCartesian 需要 3 元素数组进行索引赋值
        e0_dummy = zeros(3, 1);
        e0 = incidentFieldCartesian(uu, vv, ww, e0_dummy, Et, phr, Ep);

        % ----- 频率循环（内循环）-----
        for j = 1:N_f

            freq_j = freq_array(j);
            wave_j = c / freq_j;
            bk_j   = 2 * pi / wave_j;
            corel_j = corr / wave_j;

            % 粗糙度相关因子
            [~, cfac1_j, cfac2_j, ~, Lt_j, Nt_j] = getStandardDeviation(delstd, corel_j, wave_j);

            % 累积器（复散射场）
            sumt = 0;
            sump = 0;
            sumdt = 0;
            sumdp = 0;

            % ----- 三角面元循环 -----
            for m = 1:ntria

                % 遮挡判断（照射测试）
                ndotk = N(m, :) * R_vec;

                if iflag == 0
                    if (ilum(m) == 1 && ndotk >= 1e-5) || ilum(m) == 0

                        % 局部方向余弦
                        [u2, v2, w2, T1, T2] = directionCosines(alpha, beta, D0, m);

                        % 局部球坐标角度
                        [th2, phi2] = sphericalAngles(u2, v2, w2);

                        % 三角面元顶点的相位项（单站：因子 2*bk）
                        [Dp, Dq, Do] = phaseVerticeTriangle(x, y, z, vind, bk_j, m, u, v, w);

                        % 入射场在局部笛卡尔坐标
                        e1 = T1 * conj(e0);
                        e2_vec = T2 * e1;

                        % 入射场在局部球坐标
                        [Et2, Ep2] = incidentFieldSphericalCoordinates(th2, e2_vec, phi2);

                        % 反射系数（含频率依赖性）
                        [perp, para] = reflectionCoefficients(Rs(m), m, th2, thr, phr, ...
                            alpha(m), beta(m), freq_j, matrl);

                        % 面电流分量（局部笛卡尔坐标）
                        Jx2 = -Et2 * cos(phi2) * para + Ep2 * sin(phi2) * perp * cos(th2);
                        Jy2 = -Et2 * sin(phi2) * para - Ep2 * cos(phi2) * perp * cos(th2);

                        % 面积积分（含频率依赖的相位）
                        [DD, expDo, expDp, expDq] = areaIntegral(Dq, Dp, Do);
                        Ic = calculateIc(Dp, Dq, Do, Nt_j, Area, expDo, Co, Lt_j, DD, expDq, m, expDp);

                        % 计算散射场贡献并累加
                        [sumt, sump, sumdp, sumdt] = calculaCampos(Area, cfac2_j, corel_j, th2, wave_j, ...
                            Jy2, Ic, uu, vv, ww, phr, sumt, sump, sumdt, sumdp, m, Jx2, T1, T2);

                    end
                end
            end  % ---- 面元循环结束 ----

            % 存储复散射场的 theta 分量
            % S_complex 保存的是 coherent + diffuse 的总散射场
            S_complex(i_angle, j) = sumt;

        end  % ---- 频率循环结束 ----

    end  % ---- theta 循环结束 ----

end  % ---- phi 循环结束 ----

elapsed = toc(start_time);
fprintf('\nWideband scattering computation completed in %.2f seconds (%.1f min).\n', ...
    elapsed, elapsed/60);
fprintf('  Average: %.3f s per angle\n', elapsed / N_angles_total);

%% 7. 保存结果
fprintf('\nSaving results...\n');

% Create timestamped result directory
[resultDir, nowStr] = createResultDir('main_wideband_scattering');
resultFile = fullfile(resultDir, ['wideband_scattering_' nowStr '.mat']);

% 构建参数摘要字符串
param = sprintf(['Wideband Scattering Simulation\n', ...
    '  Model: %s\n', ...
    '  Frequency: %.2f - %.2f GHz (%d points, Δf=%.2f MHz)\n', ...
    '  Bandwidth: %.2f GHz, Range Resolution: %.4f m\n', ...
    '  Polarization: %s\n', ...
    '  Phi: %.1f:%.1f:%.1f deg (%d steps)\n', ...
    '  Theta: %.1f:%.1f:%.1f deg (%d steps)\n', ...
    '  Total angles: %d, Total frequencies: %d\n', ...
    '  Triangles: %d, Computation time: %.1f s'], ...
    inputModel, freq_start/1e9, freq_stop/1e9, N_f, freq_step/1e6, ...
    B/1e9, delta_r, pol, ...
    pstart, delp, pstop, ip, ...
    tstart, delt, tstop, it, ...
    N_angles_total, N_f, ntria, elapsed);

% 额外保存频率相关数组
wave_array = c ./ freq_array;  % 波长数组

% 保存到 .mat 文件
save(resultFile, 'S_complex', 'freq_array', 'wave_array', ...
    'phi_array', 'theta_array', 'param', 'inputModel', ...
    'c', 'B', 'delta_r', 'ip', 'it', 'N_f', 'N_angles_total', ...
    'pstart', 'pstop', 'delp', 'tstart', 'tstop', 'delt', ...
    'ntria', 'freq_start', 'freq_stop', 'pol', 'bbox_center', '-v7.3');

fprintf('  Results saved to: %s\n', resultFile);
fprintf('  Matrix size: %d angles × %d frequencies\n', N_angles_total, N_f);

%% 8. 显示摘要
fprintf('\n========================================\n');
fprintf('  Wideband Simulation Complete\n');
fprintf('========================================\n');
fprintf('  Model:           %s\n', inputModel);
fprintf('  Frequency band:  %.2f - %.2f GHz\n', freq_start/1e9, freq_stop/1e9);
fprintf('  Bandwidth:       %.2f GHz\n', B/1e9);
fprintf('  Range resolution: %.4f m\n', delta_r);
fprintf('  Frequency points: %d\n', N_f);
fprintf('  Angular range:   phi [%.1f, %.1f]°, theta [%.1f, %.1f]°\n', ...
    pstart, pstop, tstart, tstop);
fprintf('  Angle steps:     %d phi × %d theta = %d total\n', ip, it, N_angles_total);
fprintf('  Triangles:       %d\n', ntria);
fprintf('  Result file:     %s\n', resultFile);
fprintf('========================================\n');
