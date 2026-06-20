function display_po_light(stlFile, thetaInc, phiInc, thetaObs, phiObs, freqGHz)
    % DISPLAY_PO_LIGHT  Interactive PO illumination visualization on an STL model.
    %
    %   Opens a figure with:
    %     - Rotatable 3D model colored by PO scattering intensity
    %     - Control panel to change incidence/observer angles
    %     - "Compute" button for real-time update
    %     - "Save" button to export the current view as PNG
    %
    %   Colors: jet colormap (red=strong, blue=weak scatter)
    %           Black = shadowed (not illuminated)
    %
    %   Usage:
    %     >> display_po_light                           % file dialog
    %     >> display_po_light('stl_models/rocket.stl')  % with defaults
    %     >> display_po_light('stl_models/rocket.stl', 90, 0, 60, 0, 1)

    %% ---- Handle inputs ----
    if nargin < 1 || isempty(stlFile)
        [fname, pname] = uigetfile({'*.stl', 'STL Files (*.stl)'}, 'Select STL model');
        if isequal(fname, 0), disp('Cancelled.'); return; end
        stlFile = fullfile(pname, fname);
    end
    if nargin < 2, thetaInc = 90; end
    if nargin < 3, phiInc   = 0;  end
    if nargin < 4, thetaObs = 60; end
    if nargin < 5, phiObs   = 0;  end
    if nargin < 6, freqGHz  = 1;  end

    addpath('lib');

    %% ---- 1. Read STL ----
    fprintf('Reading STL: %s\n', stlFile);
    stlConverter(stlFile);
    coords = dlmread('coordinates.txt');
    facets = dlmread('facets.txt');

    % Store in a persistent data struct for callbacks
    D.r     = coords;                        % vertices [nVerts x 3]
    D.vind  = facets(:, 2:4);                % face indices [ntria x 3]
    D.ilum  = facets(:, 5);
    D.ntria = size(facets, 1);
    D.modelCenter = mean(D.r, 1);
    D.freqGHz = freqGHz;
    D.stlFile = stlFile;
    D.thetaInc = thetaInc;
    D.phiInc   = phiInc;
    D.thetaObs = thetaObs;
    D.phiObs   = phiObs;

    % Pre-compute geometry (normals validated, areas, centroids)
    fprintf('Computing geometry...\n');
    D.N = zeros(D.ntria, 3);
    D.Area = zeros(D.ntria, 1);
    D.triCentroids = zeros(D.ntria, 3);

    for i = 1:D.ntria
        A = D.r(D.vind(i,2), :) - D.r(D.vind(i,1), :);
        B = D.r(D.vind(i,3), :) - D.r(D.vind(i,2), :);
        C = D.r(D.vind(i,1), :) - D.r(D.vind(i,3), :);

        D.N(i, :) = -cross(B, A);
        Nn = norm(D.N(i, :));
        if Nn ~= 0, D.N(i, :) = D.N(i, :) / Nn; end

        d1 = norm(A); d2 = norm(B); d3 = norm(C);
        ss = 0.5 * (d1 + d2 + d3);
        D.Area(i) = sqrt(max(0, ss * (ss-d1) * (ss-d2) * (ss-d3)));

        D.triCentroids(i, :) = (D.r(D.vind(i,1), :) + D.r(D.vind(i,2), :) + D.r(D.vind(i,3), :)) / 3;
        if dot(D.N(i, :), D.triCentroids(i, :) - D.modelCenter) < 0
            D.N(i, :) = -D.N(i, :);
        end
    end

    % Pre-build face/vertex arrays for patch (geometry doesn't change)
    D.allVerts = zeros(D.ntria * 3, 3);
    D.faces = zeros(D.ntria, 3);
    for m = 1:D.ntria
        idx0 = (m-1) * 3;
        D.allVerts(idx0+1, :) = D.r(D.vind(m,1), :);
        D.allVerts(idx0+2, :) = D.r(D.vind(m,2), :);
        D.allVerts(idx0+3, :) = D.r(D.vind(m,3), :);
        D.faces(m, :) = [idx0+1, idx0+2, idx0+3];
    end

    %% ---- 2. Build figure and UI ----
    fig = figure('Name', 'PO Illumination Display', 'NumberTitle', 'off', ...
                 'Position', [50 50 1100 750], 'Color', 'w', ...
                 'CloseRequestFcn', @closeFig);

    % ---- 3D axes (left side) ----
    ax = axes('Parent', fig, 'Units', 'pixels', ...
              'Position', [30 50 800 680]);
    hold(ax, 'on');

    % Initial empty patch (will be filled on first compute)
    D.patchHandle = patch('Faces', D.faces, 'Vertices', D.allVerts, ...
                          'FaceColor', 'flat', 'EdgeColor', 'none', ...
                          'AmbientStrength', 0.3, 'DiffuseStrength', 0.8, ...
                          'SpecularStrength', 0.1, 'Parent', ax);

    % Lighting
    light('Position', [1 0 1]*10, 'Style', 'infinite', 'Color', [0.7 0.7 0.7], 'Parent', ax);
    light('Position', [-1 0 -1]*5, 'Style', 'infinite', 'Color', [0.3 0.3 0.3], 'Parent', ax);
    lighting(ax, 'gouraud');
    material(ax, 'dull');

    axis(ax, 'equal'); axis(ax, 'off');
    rotate3d(ax, 'on');
    title(ax, 'PO Illumination Map', 'FontSize', 12);

    % Store axes handle
    D.ax = ax;

    % ---- Control panel (right side) ----
    panelX = 860;  panelW = 210;
    yStart = 680;  yGap = 35;

    % Panel background
    uipanel('Parent', fig, 'Units', 'pixels', ...
            'Position', [panelX-10 20 panelW+10 yStart+30], ...
            'Title', 'Parameters', 'FontSize', 10);

    %% Input fields
    y = yStart;
    rowH = 22; labelW = 65; editW = 80; editX = panelX + labelW;

    % --- Theta Incidence ---
    uicontrol('Style', 'text', 'String', 'Theta Inc:', ...
              'Position', [panelX y labelW rowH], ...
              'HorizontalAlignment', 'right', 'FontSize', 9, ...
              'BackgroundColor', 'w');
    D.hThInc = uicontrol('Style', 'edit', 'String', num2str(thetaInc), ...
              'Position', [editX y editW rowH], ...
              'HorizontalAlignment', 'center', 'FontSize', 10, ...
              'BackgroundColor', 'w');
    uicontrol('Style', 'text', 'String', 'deg', ...
              'Position', [editX+editW+2 y 25 rowH], ...
              'HorizontalAlignment', 'left', 'FontSize', 9, ...
              'BackgroundColor', 'w');
    y = y - yGap;

    % --- Phi Incidence ---
    uicontrol('Style', 'text', 'String', 'Phi Inc:', ...
              'Position', [panelX y labelW rowH], ...
              'HorizontalAlignment', 'right', 'FontSize', 9, ...
              'BackgroundColor', 'w');
    D.hPhInc = uicontrol('Style', 'edit', 'String', num2str(phiInc), ...
              'Position', [editX y editW rowH], ...
              'HorizontalAlignment', 'center', 'FontSize', 10, ...
              'BackgroundColor', 'w');
    uicontrol('Style', 'text', 'String', 'deg', ...
              'Position', [editX+editW+2 y 25 rowH], ...
              'HorizontalAlignment', 'left', 'FontSize', 9, ...
              'BackgroundColor', 'w');
    y = y - yGap - 8;

    % --- Separator ---
    uicontrol('Style', 'text', 'String', '――――――――――', ...
              'Position', [panelX y panelW 14], ...
              'ForegroundColor', [0.6 0.6 0.6], 'FontSize', 8, ...
              'BackgroundColor', 'w');
    y = y - 22;

    % --- Theta Observer ---
    uicontrol('Style', 'text', 'String', 'Theta Obs:', ...
              'Position', [panelX y labelW rowH], ...
              'HorizontalAlignment', 'right', 'FontSize', 9, ...
              'BackgroundColor', 'w');
    D.hThObs = uicontrol('Style', 'edit', 'String', num2str(thetaObs), ...
              'Position', [editX y editW rowH], ...
              'HorizontalAlignment', 'center', 'FontSize', 10, ...
              'BackgroundColor', 'w');
    uicontrol('Style', 'text', 'String', 'deg', ...
              'Position', [editX+editW+2 y 25 rowH], ...
              'HorizontalAlignment', 'left', 'FontSize', 9, ...
              'BackgroundColor', 'w');
    y = y - yGap;

    % --- Phi Observer ---
    uicontrol('Style', 'text', 'String', 'Phi Obs:', ...
              'Position', [panelX y labelW rowH], ...
              'HorizontalAlignment', 'right', 'FontSize', 9, ...
              'BackgroundColor', 'w');
    D.hPhObs = uicontrol('Style', 'edit', 'String', num2str(phiObs), ...
              'Position', [editX y editW rowH], ...
              'HorizontalAlignment', 'center', 'FontSize', 10, ...
              'BackgroundColor', 'w');
    uicontrol('Style', 'text', 'String', 'deg', ...
              'Position', [editX+editW+2 y 25 rowH], ...
              'HorizontalAlignment', 'left', 'FontSize', 9, ...
              'BackgroundColor', 'w');
    y = y - yGap - 10;

    % --- Compute button ---
    D.hCompute = uicontrol('Style', 'pushbutton', 'String', 'Compute', ...
              'Position', [panelX y panelW 32], ...
              'FontSize', 11, 'FontWeight', 'bold', ...
              'BackgroundColor', [0.3 0.6 1.0], 'ForegroundColor', 'w', ...
              'Callback', @(src, evt) computeAndUpdate());
    y = y - 42;

    % --- Save button ---
    D.hSave = uicontrol('Style', 'pushbutton', 'String', 'Save Image', ...
              'Position', [panelX y panelW 28], ...
              'FontSize', 10, ...
              'BackgroundColor', [0.5 0.8 0.5], ...
              'Callback', @(src, evt) saveCurrentView());
    y = y - 40;

    % --- Info text ---
    D.hInfo = uicontrol('Style', 'text', ...
              'Position', [panelX 20 panelW 60], ...
              'String', 'Ready.', ...
              'HorizontalAlignment', 'left', 'FontSize', 8, ...
              'BackgroundColor', 'w', 'ForegroundColor', [0.4 0.4 0.4]);

    %% Store data and do initial compute
    guidata(fig, D);
    computeAndUpdate();

    %% ==================== Nested Callback Functions ====================

    function computeAndUpdate()
        % Read current parameters from UI
        D = guidata(fig);

        thetaInc = str2double(get(D.hThInc, 'String'));
        phiInc   = str2double(get(D.hPhInc, 'String'));
        thetaObs = str2double(get(D.hThObs, 'String'));
        phiObs   = str2double(get(D.hPhObs, 'String'));

        if isnan(thetaInc) || isnan(phiInc) || isnan(thetaObs) || isnan(phiObs)
            set(D.hInfo, 'String', 'Error: invalid angle value');
            return;
        end

        set(D.hInfo, 'String', 'Computing...'); drawnow;

        % ---- Compute illumination ----
        rad = pi / 180;
        sti = sin(thetaInc*rad); cti = cos(thetaInc*rad);
        spi = sin(phiInc*rad);   cpi = cos(phiInc*rad);
        D_inc = [sti*cpi; sti*spi; cti];

        sto = sin(thetaObs*rad); cto = cos(thetaObs*rad);
        spo = sin(phiObs*rad);   cpo = cos(phiObs*rad);
        D_obs = [sto*cpo; sto*spo; cto];

        ndotk_inc = D.N * D_inc;
        lit_std = (D.ilum == 1 & ndotk_inc >= 1e-5) | (D.ilum == 0);
        depthDiff = (D.triCentroids - D.modelCenter) * D_inc;
        lit_final = lit_std & (depthDiff >= 0 | ndotk_inc >= 0.01);

        % ---- Compute scattering intensity ----
        Et = 1; Ep = 0;
        uui = cti*cpi; vvi = cti*spi; wwi = -sti;
        E_inc = [uui*Et - spi*Ep; vvi*Et + cpi*Ep; wwi*Et];
        k_i = D_inc;
        H_inc = cross(k_i, E_inc);
        H_inc = H_inc / (norm(H_inc) + 1e-30);

        intensity = zeros(D.ntria, 1);
        for m = 1:D.ntria
            if ~lit_final(m), continue; end
            Js = 2 * cross(D.N(m, :), H_inc);
            Js_mag = norm(Js);
            if Js_mag < 1e-30, continue; end
            Js_proj = Js - dot(Js, D_obs) * D_obs;
            scat_mag = norm(Js_proj);
            localFactor = max(0, ndotk_inc(m));
            intensity(m) = D.Area(m) * scat_mag * localFactor;
        end

        maxInt = max(intensity);
        if maxInt > 0
            intensity_norm = intensity / maxInt;
        else
            intensity_norm = intensity;
        end

        % ---- Update patch colors ----
        faceColors = zeros(D.ntria, 3);
        for m = 1:D.ntria
            if lit_final(m) && intensity_norm(m) > 1e-6
                faceColors(m, :) = jetColor(intensity_norm(m));
            else
                faceColors(m, :) = [0 0 0];
            end
        end

        set(D.patchHandle, 'FaceVertexCData', faceColors);

        % ---- Update colormap ----
        colormap(D.ax, jet(256));
        cmap = colormap(D.ax);
        cmap(1, :) = [0 0 0];
        colormap(D.ax, cmap);
        clim(D.ax, [0 1]);

        % Update colorbar (remove old, add new)
        delete(findobj(fig, 'Type', 'ColorBar'));
        cbar = colorbar(D.ax);
        cbar.Label.String = 'Normalized Scattering Intensity';
        cbar.Label.FontSize = 9;

        % ---- Set camera to observer perspective ----
        obsDist = max(max(D.r) - min(D.r)) * 3;
        campos(D.ax, D_obs' * obsDist);
        camtarget(D.ax, D.modelCenter);
        camup(D.ax, [0 0 1]);
        camproj(D.ax, 'perspective');
        camva(D.ax, 15);

        % ---- Update title ----
        title(D.ax, sprintf(['PO Illumination  |  Drag to rotate\n' ...
               'Inc: \\theta=%.0f°, \\phi=%.0f°  →  ' ...
               'Obs: \\theta=%.0f°, \\phi=%.0f°'], ...
               thetaInc, phiInc, thetaObs, phiObs), 'FontSize', 12);

        % ---- Update info ----
        set(D.hInfo, 'String', sprintf('Lit: %d/%d (%.1f%%)', ...
                sum(lit_final), D.ntria, 100*sum(lit_final)/D.ntria));

        % Store updated values
        D.thetaInc = thetaInc; D.phiInc = phiInc;
        D.thetaObs = thetaObs; D.phiObs = phiObs;
        guidata(fig, D);
    end

    function saveCurrentView()
        D = guidata(fig);
        [~, modelName] = fileparts(D.stlFile);
        outName = sprintf('po_light_%s_incTH%d_PH%d_obsTH%d_PH%d.png', ...
                          modelName, round(D.thetaInc), round(D.phiInc), ...
                          round(D.thetaObs), round(D.phiObs));
        saveas(fig, outName);
        set(D.hInfo, 'String', sprintf('Saved: %s', outName));
        fprintf('Figure saved: %s\n', outName);
    end

    function closeFig(~, ~)
        delete(fig);
    end
end

%% ==================== Helper: jet colormap interpolation ====================
function rgb = jetColor(t)
    t = max(0, min(1, t));
    r = interp1([0 0.25 0.5 0.75 1], [0 0 1 1 0.5], t, 'linear');
    g = interp1([0 0.25 0.5 0.75 1], [0 1 1 0 0],   t, 'linear');
    b = interp1([0 0.25 0.5 0.75 1], [0.5 1 0 0 0], t, 'linear');
    rgb = [r g b];
end
