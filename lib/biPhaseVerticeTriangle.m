function [Dp, Dq, Do] = biPhaseVerticeTriangle(x, y, z, vind, bk, m, u, v, w, ui, vi, wi)
    % BIPHASEVERTICETRIANGLE Compute phase terms at triangle vertices for bistatic RCS.
    %
    %   Bistatic: uses factor bk with sum of incidence and observation directions.
    %
    %   Input:
    %       x, y, z     - vertex coordinate vectors
    %       vind        - vertex index table
    %       bk          - wave number
    %       m           - triangle index
    %       u, v, w     - observation direction cosines
    %       ui, vi, wi  - incidence direction cosines
    %   Output:
    %       Dp, Dq, Do - phase terms at triangle vertices

    Dp = bk * ((x(vind(m, 1)) - x(vind(m, 3))) * (u + ui) + ...
               (y(vind(m, 1)) - y(vind(m, 3))) * (v + vi) + ...
               (z(vind(m, 1)) - z(vind(m, 3))) * (w + wi));

    Dq = bk * ((x(vind(m, 2)) - x(vind(m, 3))) * (u + ui) + ...
               (y(vind(m, 2)) - y(vind(m, 3))) * (v + vi) + ...
               (z(vind(m, 2)) - z(vind(m, 3))) * (w + wi));

    Do = bk * (x(vind(m, 3)) * (u + ui) + y(vind(m, 3)) * (v + vi) + z(vind(m, 3)) * (w + wi));
end
