function plotName = finalPlot(ip, it, phi, wave, theta, Lmin, Lmax, Sth, Sph, U, V, nowStr, inputModel, mode, resultDir)
    % FINALPLOT Generate RCS plots: theta-cut, phi-cut, or contour.
    %
    %   Input:
    %       ip, it     - number of phi and theta steps
    %       phi, theta - angle arrays (ip x it)
    %       wave       - wavelength (m)
    %       Lmin, Lmax - axis limits
    %       Sth, Sph   - RCS arrays (ip x it)
    %       U, V       - direction cosine arrays
    %       nowStr     - timestamp string
    %       inputModel - model name
    %       mode       - 'Monostatic' or 'Bistatic'
    %       resultDir  - (optional) result subdirectory path; defaults to 'results'
    %   Output:
    %       plotName - path to saved plot

    if nargin < 15 || isempty(resultDir)
        resultDir = 'results';
        if ~exist(resultDir, 'dir'), mkdir(resultDir); end
    end

    fig = figure('Visible', 'on');

    if ip == 1
        % Single phi value: plot vs theta
        sgtitle(sprintf('RCS Simulation IR Signature - %s', mode));
        title(sprintf('target: %s   solid: theta     dashed: phi     phi= %.1f    wave (m): %.6f', ...
              inputModel, phi(1, 1), round(wave, 6)));
        xlabel(sprintf('%s Angle, theta (deg)', mode));
        ylabel('RCS (dBsm)');
        axis([min(theta(:)), max(theta(:)), Lmin, Lmax]);
        hold on;
        plot(theta(1, :), Sth(1, :), 'b-', 'LineWidth', 1.5);
        plot(theta(1, :), Sph(1, :), 'r--', 'LineWidth', 1.5);
        legend('S_{theta}', 'S_{phi}', 'Location', 'best');
        grid on;

    elseif it == 1
        % Single theta value: plot vs phi
        sgtitle(sprintf('RCS Simulation IR Signature - %s', mode));
        title(sprintf('target: %s   solid: theta     dashed: phi     theta= %.1f    wave (m): %.6f', ...
              inputModel, theta(1, 1), round(wave, 6)));
        xlabel(sprintf('%s Angle, phi (deg)', mode));
        ylabel('RCS (dBsm)');
        axis([min(phi(:)), max(phi(:)), Lmin, Lmax]);
        hold on;
        plot(phi(:, 1), Sth(:, 1), 'b-', 'LineWidth', 1.5);
        plot(phi(:, 1), Sph(:, 1), 'r--', 'LineWidth', 1.5);
        legend('S_{theta}', 'S_{phi}', 'Location', 'best');
        grid on;

    elseif ip > 1 && it > 1
        % 2D contour plot
        sgtitle(sprintf('RCS Simulation IR Signature - %s', mode));

        % RCS-theta contour
        subplot(2, 3, 2);
        if strcmp(mode, 'Monostatic')
            contour(U, V, Sth);
        else
            contour(U, V, Sth, [-20, 0]);
        end
        title('RCS-theta');
        xlabel('U');
        ylabel('V');
        axis square;
        colorbar;

        % RCS-phi contour
        subplot(2, 3, 5);
        if strcmp(mode, 'Monostatic')
            contour(U, V, Sph);
        else
            contour(U, V, Sph, [-20, 0]);
        end
        title('RCS-phi');
        xlabel('U');
        ylabel('V');
        axis square;
        colorbar;
    end

    % Save plot
    plotName = fullfile(resultDir, ['temp_' nowStr '.png']);
    saveas(fig, plotName);
end
