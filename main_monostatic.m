% MAIN_MONOSTATIC  Monostatic Radar Cross Section (RCS) Simulation
%
%   Physical Optics (PO) based RCS computation for arbitrary 3D targets
%   defined by STL triangular mesh models.
%
%   This script:
%     1. Reads simulation parameters from input_files/input_data_file_monostatic.dat
%     2. Converts the specified STL model to coordinates and facets
%     3. Computes monostatic RCS using the Physical Optics method
%     4. Generates plots and saves results to the results/ directory
%
%   Usage:
%     >> main_monostatic

clear; clc;

% Add library path
addpath('lib');

fprintf('========================================\n');
fprintf('  Open-RCS: Monostatic RCS Simulation\n');
fprintf('========================================\n\n');

%% 1. Load simulation parameters
fprintf('Loading parameters...\n');
paramList = getParamsFromFile('monostatic');

inputModel = paramList{1};
freq = paramList{2};        % Hz
corr = paramList{3};        % m
delstd = paramList{4};      % m
ipol = paramList{5};        % 0=TM-z, 1=TE-z
rs = paramList{6};          % 0=PEC, 1=material
pstart = paramList{7};      % deg
pstop = paramList{8};       % deg
delp = paramList{9};        % deg
tstart = paramList{10};     % deg
tstop = paramList{11};      % deg
delt = paramList{12};       % deg
matrlpath = paramList{end}; % material file path

MATERIALESPECIFICO = 1;

%% 2. Extract geometry data
fprintf('Processing geometry...\n');
coordinatesData = extractCoordinatesData(rs);

x = coordinatesData{1};
y = coordinatesData{2};
z = coordinatesData{3};
xpts = coordinatesData{4};
ypts = coordinatesData{5};
zpts = coordinatesData{6};
nverts = coordinatesData{7};
nfc = coordinatesData{8};
node1 = coordinatesData{9};
node2 = coordinatesData{10};
node3 = coordinatesData{11};
iflag = coordinatesData{12};
ilum = coordinatesData{13};
Rs = coordinatesData{14};
ntria = coordinatesData{15};
vind = coordinatesData{16};
r = coordinatesData{17};

%% 3. Load material properties if needed
matrl = {};
if rs == MATERIALESPECIFICO
    try
        matrl = getEntrysFromMatrlFile(ntria, matrlpath);
        fprintf('Material file loaded: %s\n', matrlpath);
    catch ME
        error('Material file error: %s', ME.message);
    end
end

%% 4. Physical parameters
wave = 3e8 / freq;
corel = corr / wave;
[bk, cfac1, cfac2, rad, Lt, Nt] = getStandardDeviation(delstd, corel, wave);
[pol, Et, Ep] = getPolarization(ipol);
Co = 1;  % wave amplitude at all vertices

fprintf('Frequency: %.2f GHz\n', freq/1e9);
fprintf('Wavelength: %.6f m\n', wave);
fprintf('Polarization: %s\n', pol);
fprintf('Triangles: %d\n', ntria);

%% 5. Pre-compute geometry arrays
fprintf('Computing triangle geometry...\n');
[Area, alpha, beta, N, d, ip, it] = calculateValues(pstart, pstop, delp, tstart, tstop, delt, ntria, rad);
[N, d, Area, beta, alpha] = productVector(ntria, N, r, d, Area, alpha, beta, vind);
[phi, theta, U, V, W, e0, Sth, Sph] = otherVectorComponents(ip, it);

% Pre-compute model center and triangle centroids for self-shadowing check
modelCenter = mean(r, 1);
triCentroids = zeros(ntria, 3);
for m = 1:ntria
    triCentroids(m, :) = (r(vind(m, 1), :) + r(vind(m, 2), :) + r(vind(m, 3), :)) / 3;
end

fprintf('Angular grid: %d phi x %d theta points\n', ip, it);

%% 6. Main RCS computation loop
fprintf('Computing RCS...\n');
tic;

for i1 = 1:ip
    for i2 = 1:it
        % Observation angles
        phi(i1, i2) = pstart + (i1 - 1) * delp;
        phr = phi(i1, i2) * rad;
        theta(i1, i2) = tstart + (i2 - 1) * delt;
        thr = theta(i1, i2) * rad;

        % Global direction cosines
        [U, V, W, D0, uu, vv, ww, u, v, w] = globalAngles(U, V, W, thr, phr, i1, i2);

        % Radial unit vector
        R_vec = [u; v; w];

        % Incident field in global Cartesian coordinates
        e0 = incidentFieldCartesian(uu, vv, ww, e0, Et, phr, Ep);

        % Accumulators
        sumt = 0;
        sump = 0;
        sumdt = 0;
        sumdp = 0;

        % Loop over all triangles
        for m = 1:ntria
            % Illumination test
            ndotk = N(m, :) * R_vec;

            if iflag == 0
                if (ilum(m) == 1 && ndotk >= 1e-5) || ilum(m) == 0
                    % Self-shadowing check: triangles on the back face of the
                    % model (relative to illumination direction) with marginal
                    % N·D > 0 may be incorrectly lit due to body curvature.
                    % Require either: (a) triangle is on the front face, or
                    % (b) normal is strongly forward-facing.
                    depthDiff = dot(triCentroids(m, :) - modelCenter, R_vec);
                    if depthDiff >= 0 || ndotk >= 0.01
                    % Local direction cosines
                    [u2, v2, w2, T1, T2] = directionCosines(alpha, beta, D0, m);

                    % Local spherical angles
                    [th2, phi2] = sphericalAngles(u2, v2, w2);

                    % Phase at triangle vertices (monostatic: factor 2*bk)
                    [Dp, Dq, Do] = phaseVerticeTriangle(x, y, z, vind, bk, m, u, v, w);

                    % Incident field in local Cartesian coordinates
                    e1 = T1 * conj(e0);
                    e2_vec = T2 * e1;

                    % Incident field in local spherical coordinates
                    [Et2, Ep2] = incidentFieldSphericalCoordinates(th2, e2_vec, phi2);

                    % Reflection coefficients
                    [perp, para] = reflectionCoefficients(Rs(m), m, th2, thr, phr, alpha(m), beta(m), freq, matrl);

                    % Surface current components in local Cartesian
                    Jx2 = -Et2 * cos(phi2) * para + Ep2 * sin(phi2) * perp * cos(th2);
                    Jy2 = -Et2 * sin(phi2) * para - Ep2 * cos(phi2) * perp * cos(th2);

                    % Area integral
                    [DD, expDo, expDp, expDq] = areaIntegral(Dq, Dp, Do);
                    Ic = calculateIc(Dp, Dq, Do, Nt, Area, expDo, Co, Lt, DD, expDq, m, expDp);

                    % Compute scattered fields
                    [sumt, sump, sumdp, sumdt] = calculaCampos(Area, cfac2, corel, th2, wave, ...
                        Jy2, Ic, uu, vv, ww, phr, sumt, sump, sumdt, sumdp, m, Jx2, T1, T2);
                    end  % self-shadowing check
                end
            end
        end

        % Compute RCS for this angular position
        [Sth, Sph] = calculateSthSph(cfac1, sumt, sump, sumdt, wave, Sth, Sph, i1, i2, sumdp);
    end

    % Progress indicator
    if mod(i1, max(1, floor(ip/10))) == 0
        fprintf('  Progress: %d/%d phi angles (%.1f%%)\n', i1, ip, 100*i1/ip);
    end
end

elapsed = toc;
fprintf('RCS computation completed in %.2f seconds.\n', elapsed);

%% 7. Compute axis limits and generate outputs
SthPlot = Sth;
SphPlot = Sph;
[Lmax, Lmin] = parametrosGrafico(SthPlot, SphPlot);
Lmax=30; Lmin=-40;
fprintf('RCS range: %.1f to %.1f dBsm\n', Lmin, Lmax);

%% 8. Generate result files
fprintf('Generating results...\n');

% Create timestamped result directory
[resultDir, nowStr] = createResultDir('main_monostatic');

% Triangle model figure
setFontOption();
% figName = plotTriangleModel(inputModel, vind, x, y, z, xpts, ypts, zpts, nverts, ntria, node1, node2, node3, nfc, resultDir);

% Parameter summary
param = plotParameters('Monostatic', freq, wave, corr, delstd, pol, ntria, pstart, pstop, delp, tstart, tstop, delt);

% Data file
[nowStr, fileName] = generateResultFiles(theta, Sth, phi, Sph, param, ip, resultDir);

% RCS plot
plotName = finalPlot(ip, it, phi, wave, theta, Lmin, Lmax, SthPlot, SphPlot, U, V, nowStr, inputModel, 'Monostatic', resultDir);

%% 9. Display summary
fprintf('\n========================================\n');
fprintf('  Simulation Complete\n');
fprintf('========================================\n');
fprintf('  Model:    %s\n', inputModel);
fprintf('  Plot:     %s\n', plotName);
fprintf('  Figure:   %s\n', figName);
fprintf('  Data:     %s\n', fileName);
fprintf('========================================\n');
