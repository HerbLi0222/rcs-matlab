function [U, V, W, D0, uu, vv, ww, u, v, w] = globalAngles(U, V, W, thr, phr, i1, i2)
    % GLOBALANGLES Compute global direction cosines and angular components.
    %
    %   Input:
    %       U, V, W - direction cosine arrays (updated in-place)
    %       thr     - theta angle (radians)
    %       phr     - phi angle (radians)
    %       i1, i2  - array indices
    %
    %   Output:
    %       U, V, W  - direction cosine arrays
    %       D0       - direction cosine vector [u; v; w]
    %       uu, vv, ww - theta-derivative direction cosines
    %       u, v, w  - direction cosines (scalars)

    u = sin(thr) * cos(phr);
    v = sin(thr) * sin(phr);
    w = cos(thr);

    U(i1, i2) = u;
    V(i1, i2) = v;
    W(i1, i2) = w;

    D0 = [u; v; w];

    uu = cos(thr) * cos(phr);
    vv = cos(thr) * sin(phr);
    ww = -sin(thr);
end
