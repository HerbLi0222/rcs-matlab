function [Lmax, Lmin] = parametrosGrafico(Sth, Sph)
    % PARAMETROSGRAFICO Calculate plot axis limits from RCS data.
    %
    %   Input:
    %       Sth, Sph - RCS arrays (dBsm)
    %   Output:
    %       Lmax - upper axis limit (next multiple of 5 above max)
    %       Lmin - lower axis limit (minimum of all values)

    Smax = max(max(Sth(:)), max(Sph(:)));
    Lmax = (floor(Smax / 5) + 1) * 5;
    Lmin = min(min(Sth(:)), min(Sph(:)));
end
