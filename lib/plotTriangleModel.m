function figName = plotTriangleModel(inputModel, vind, x, y, z, xpts, ypts, zpts, nverts, ntria, node1, node2, node3, nfc, resultDir)
    % PLOTTRIANGLEMODEL Plot 3D wireframe of the triangular mesh model.
    %
    %   Input:
    %       inputModel - model filename
    %       vind       - vertex index table
    %       x, y, z    - vertex coordinate vectors
    %       xpts, ypts, zpts - vertex coordinates
    %       nverts     - number of vertices
    %       ntria      - number of triangles
    %       node1, node2, node3 - facet node indices
    %       nfc        - facet numbers
    %       resultDir  - (optional) result subdirectory path; defaults to 'results'
    %   Output:
    %       figName - path to saved figure

    if nargin < 15 || isempty(resultDir)
        resultDir = 'results';
        if ~exist(resultDir, 'dir'), mkdir(resultDir); end
    end

    fig = figure('Visible', 'on');
    ax = axes('Parent', fig);

    hold(ax, 'on');

    % Plot each triangle
    for i = 1:ntria
        X = [x(vind(i, 1)), x(vind(i, 2)), x(vind(i, 3)), x(vind(i, 1))];
        Y = [y(vind(i, 1)), y(vind(i, 2)), y(vind(i, 3)), y(vind(i, 1))];
        Z = [z(vind(i, 1)), z(vind(i, 2)), z(vind(i, 3)), z(vind(i, 1))];
        plot3(ax, X, Y, Z, 'b-', 'LineWidth', 0.5);
    end

    xlabel(ax, 'x');
    ylabel(ax, 'y');
    zlabel(ax, 'z');

    title(ax, sprintf('Triangle Model of Target: %s', inputModel));

    % Equal aspect ratio
    xmax = max(xpts); xmin = min(xpts);
    ymax = max(ypts); ymin = min(ypts);
    zmax = max(zpts); zmin = min(zpts);

    x_range = xmax - xmin;
    y_range = ymax - ymin;
    z_range = zmax - zmin;
    max_range = max([x_range, y_range, z_range]);

    xlim(ax, [xmin, xmin + max_range]);
    ylim(ax, [ymin, ymin + max_range]);
    zlim(ax, [zmin, zmin + max_range]);

    axis(ax, 'equal');
    view(ax, 3);
    grid(ax, 'on');

    % Save figure
    nowStr = datestr(now, 'yyyymmddHHMMSS');
    figName = fullfile(resultDir, ['temp_' nowStr '.jpg']);

    saveas(fig, figName);
end
