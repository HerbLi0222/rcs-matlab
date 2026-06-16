function g = GFunc(n, w)
    % GFUNC Recursive G function for Physical Optics phase integration.
    %   G(n, w) = (exp(j*w) - 1) / (j*w)  for n = 0
    %   G(n, w) = (exp(j*w) - n * G(n-1, w)) / (j*w)  for n > 0
    %
    %   Input:
    %       n - order (non-negative integer)
    %       w - argument (real scalar)
    %   Output:
    %       g - complex result

    jw = 1i * w;
    g = (exp(jw) - 1) / jw;

    if n > 0
        for m = 1:n
            go = g;
            g = (exp(jw) - n * go) / jw;
        end
    end
end
