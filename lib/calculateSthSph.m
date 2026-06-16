function [Sth, Sph] = calculateSthSph(cfac1, sumt, sump, sumdt, wave, Sth, Sph, i1, i2, sumdp)
    % CALCULATESTHSPH Compute final RCS values in dBsm.
    %
    %   RCS = 10 * log10(4*pi * [cfac1*|sum|^2 + sqrt(1-cfac1^2)*sumDiff] / lambda^2)
    %
    %   A small epsilon (1e-10) is added inside the log to avoid log(0).
    %
    %   Input:
    %       cfac1 - coherent scattering factor
    %       sumt, sump - coherent field sums (theta, phi)
    %       sumdt, sumdp - diffuse field sums
    %       wave  - wavelength (m)
    %       Sth, Sph - RCS output arrays (updated in place)
    %       i1, i2 - array indices

    Sth(i1, i2) = 10 * log10(4 * pi * cfac1 * (abs(sumt)^2 + sqrt(1 - cfac1^2) * sumdt) / wave^2 + 1e-10);
    Sph(i1, i2) = 10 * log10(4 * pi * cfac1 * (abs(sump)^2 + sqrt(1 - cfac1^2) * sumdp) / wave^2 + 1e-10);
end
