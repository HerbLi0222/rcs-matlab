function [Dp, Dq, Do] = phaseVerticeTriangle(x, y, z, vind, bk, m, u, v, w)
    % PHASEVERTICETRIANGLE Compute phase terms at triangle vertices for monostatic RCS.
    %
    %   Monostatic: uses factor 2*bk (round-trip phase).
    %
    %   Input:
    %       x, y, z - vertex coordinate vectors
    %       vind    - vertex index table
    %       bk      - wave number (2*pi/wavelength)
    %       m       - triangle index
    %       u, v, w - direction cosines
    %   Output:
    %       Dp, Dq, Do - phase terms at triangle vertices

    Dp = 2 * bk * ((x(vind(m, 1)) - x(vind(m, 3))) * u + ...
                   (y(vind(m, 1)) - y(vind(m, 3))) * v + ...
                   (z(vind(m, 1)) - z(vind(m, 3))) * w);

    Dq = 2 * bk * ((x(vind(m, 2)) - x(vind(m, 3))) * u + ...
                   (y(vind(m, 2)) - y(vind(m, 3))) * v + ...
                   (z(vind(m, 2)) - z(vind(m, 3))) * w);

    Do = 2 * bk * (x(vind(m, 3)) * u + y(vind(m, 3)) * v + z(vind(m, 3)) * w);
end
