function y = clamp(x, a, b)
    y = reshape(max(a, min(b, x(:))), size(x));
end