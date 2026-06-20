% MAIN_DISPLAY_PO_LIGHT  Main script for PO illumination visualization
%
%   Configure parameters below and run to visualize which parts of the
%   STL model are illuminated by a plane wave and how strongly they
%   scatter toward the observer.
%
%   Illuminated triangles: jet colormap (red=strong, blue=weak scatter)
%   Shadowed triangles:   black
%   Viewpoint:            observer's perspective
%
%   Usage:
%     >> main_display_po_light

clear; clc;

fprintf('========================================\n');
fprintf('  PO Illumination Visualization\n');
fprintf('========================================\n\n');

%% ==================== CONFIGURATION ====================

% --- Model ---
stlFile = 'stl_models/rocket.stl';   % STL model path

% --- Radar Frequency (for display info only) ---
freqGHz = 1;                          % GHz

% --- Incident wave direction (where the radar is) ---
%     theta = 0: +z direction (nose-on for rocket)
%     theta = 90: broadside
%     theta = 180: -z direction (tail-on)
thetaInc = 90;                         % deg, 0~180
phiInc   = 0;                         % deg, 0~360

% --- Observer direction (where you look from) ---
%     Single view: set thetaObs / phiObs
%     Multi-view:  set thetaObsArray / phiObsArray (generates one figure each)
thetaObs = 60;                        % deg, 0~180
phiObs   = 0;                        % deg, 0~360

% --- Multi-view mode (set to true to sweep observer angles) ---
multiView = false;
% Observer sweep parameters (only used if multiView = true)
thetaObsArray = 0:30:180;            % sweep theta
phiObsFixed   = 0;                    % fixed phi for sweep

%% ==================== RUN ====================

addpath('lib');

if ~multiView
    %% Single view
    fprintf('Incident:  theta = %.0f deg, phi = %.0f deg\n', thetaInc, phiInc);
    fprintf('Observer:  theta = %.0f deg, phi = %.0f deg\n', thetaObs, phiObs);
    fprintf('Model:     %s\n\n', stlFile);

    display_po_light(stlFile, thetaInc, phiInc, thetaObs, phiObs, freqGHz);
else
    %% Multi-view sweep
    nViews = length(thetaObsArray);
    fprintf('Incident:  theta = %.0f deg, phi = %.0f deg\n', thetaInc, phiInc);
    fprintf('Observer sweep: theta = [%s] deg, phi = %.0f deg\n', ...
            num2str(thetaObsArray), phiObsFixed);
    fprintf('Model:     %s\n', stlFile);
    fprintf('Generating %d views...\n\n', nViews);

    for v = 1:nViews
        thObs = thetaObsArray(v);
        fprintf('[%d/%d] Observer: theta = %.0f deg, phi = %.0f deg\n', ...
                v, nViews, thObs, phiObsFixed);
        display_po_light(stlFile, thetaInc, phiInc, thObs, phiObsFixed, freqGHz);
    end
end

fprintf('\n========================================\n');
fprintf('  Done.\n');
fprintf('========================================\n');
