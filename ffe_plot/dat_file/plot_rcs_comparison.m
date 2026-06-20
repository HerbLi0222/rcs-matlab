% plot_rcs_comparison.m
% Compare monostatic RCS from two sources in the same figure:
%   1. FEKO simulation (result_rocket.ffe) — POSTFEKO bistatic far-field export
%   2. Custom simulator  (temp_*.dat)      — TM-z polarization monostatic
%
% Both data sets are at 1 GHz, theta = 0°–180°, phi = 0° (monostatic).
% The .ffe file contains bistatic data; the monostatic RCS is extracted
% along the theta_obs = theta_inc diagonal.

clear; clc; close all;

%% ----- Path setup -----
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
addpath(fullfile(script_dir, '..'));

%% ----- Locate input files -----
ffe_file = fullfile(script_dir, 'result_rocket_po.ffe');
dat_file = dir(fullfile(script_dir, 'temp_*.dat'));
assert(exist(ffe_file, 'file'), 'FFE file not found: %s', ffe_file);
assert(~isempty(dat_file), 'No temp_*.dat file found in: %s', script_dir);
dat_file = fullfile(script_dir, dat_file(1).name);

%% ----- 1. Load FEKO .ffe data & extract monostatic diagonal -----
ffe = ffe_load_data(ffe_file);

theta_obs = ffe.theta(1, :)';      % [181×1], -180° to 180°
theta_inc = ffe.theta_inc(:);      % [181×1], 0° to 180°
RCS_total_db = 10 * log10(max(ffe.RCS_total, 1e-30));

mono_theta = zeros(length(theta_inc), 1);
mono_rcs   = zeros(length(theta_inc), 1);
for i = 1:length(theta_inc)
    [~, obs_idx] = min(abs(theta_obs - theta_inc(i)));
    mono_theta(i) = theta_inc(i);
    mono_rcs(i)   = RCS_total_db(obs_idx, i);
end

%% ----- 2. Load simulator .dat data -----
sim = ffe_load_data(dat_file);
sim_theta = sim.theta(1, :)';
sim_rcs   = sim.Sth(1, :)';

%% ----- 3. Plot -----
figure('Name', 'Monostatic RCS Comparison', 'NumberTitle', 'off', ...
       'Position', [100, 100, 900, 550]);
hold on;
plot(mono_theta, mono_rcs, 'b-', 'LineWidth', 1, ...
     'DisplayName', sprintf('FEKO %s', ffe.fileName));
plot(sim_theta, sim_rcs, 'r-', 'LineWidth', 1, ...
     'DisplayName', sprintf('Simulator %s', sim.fileName));
hold off;

xlabel('\theta (deg)', 'FontSize', 12);
ylabel('RCS (dBsm)', 'FontSize', 12);
title(sprintf('Monostatic RCS Comparison @ %.2f GHz', ffe.freq_ghz), 'FontSize', 13);
legend('Location', 'best', 'FontSize', 10);
grid on; box on;
xlim([0 180]);

fprintf('Done.\n');
