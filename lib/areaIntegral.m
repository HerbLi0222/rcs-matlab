function [DD, expDo, expDp, expDq] = areaIntegral(Dq, Dp, Do)
    % AREAINTEGRAL Compute pre-factors for the area integral.
    %
    %   Input:
    %       Dp, Dq, Do - phase terms at triangle vertices
    %   Output:
    %       DD    - Dq - Dp
    %       expDo - exp(1i * Do)
    %       expDp - exp(1i * Dp)
    %       expDq - exp(1i * Dq)

    DD = Dq - Dp;
    expDo = exp(1i * Do);
    expDp = exp(1i * Dp);
    expDq = exp(1i * Dq);
end
