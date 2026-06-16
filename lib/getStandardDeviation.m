function [bk, cfac1, cfac2, rad, Lt, Nt] = getStandardDeviation(delstd, corel, wave)
    % GETSTANDARDDEVIATION Compute roughness-related parameters for RCS calculation.
    %
    %   Input:
    %       delstd - standard deviation of surface roughness (m)
    %       corel  - correlation distance normalized to wavelength
    %       wave   - wavelength (m)
    %   Output:
    %       bk    - wave number (2*pi/wave)
    %       cfac1 - coherent factor: exp(-4*bk^2*delsq)
    %       cfac2 - diffuse factor: 4*pi*(bk*corel)^2*delsq
    %       rad   - degrees to radians conversion (pi/180)
    %       Lt    - Taylor series threshold
    %       Nt    - number of Taylor series terms

    delsq = delstd^2;
    bk = 2 * pi / wave;
    cfac1 = exp(-4 * bk^2 * delsq);
    cfac2 = 4 * pi * (bk * corel)^2 * delsq;
    rad = pi / 180;
    Lt = 1e-5;   % Taylor series region
    Nt = 5;       % number of terms in Taylor series
end
