function Es_corrected = alignPhaseReference(Es_f, freq_array, theta_deg, phi_deg, r0, c)
% ALIGNPHASEREFERENCE  Phase-ramp range alignment for HRRP.
%
%   Shifts the phase reference from the coordinate origin to point r0
%   by applying a frequency-domain phase ramp:
%       E'(f) = E(f) * exp(-j * 4*pi*f/c * r0 . k_hat)
%
%   This is the standard radar range-alignment technique. It does NOT
%   require modifying geometry coordinates or re-running scattering.
%
%   Input:
%       Es_f       - complex scattered field, 1 x N_f
%       freq_array - frequency samples (Hz), 1 x N_f
%       theta_deg  - observation theta angle (deg)
%       phi_deg    - observation phi angle (deg)
%       r0         - new reference point [x0; y0; z0] (m), 3 x 1
%       c          - speed of light (m/s)
%   Output:
%       Es_corrected - phase-aligned complex field, 1 x N_f

if nargin < 6, c = 3e8; end
if nargin < 5 || isempty(r0), r0 = [0; 0; 0]; end

% Skip if reference point is at origin
if norm(r0) < 1e-10
    Es_corrected = Es_f;
    return;
end

% Radar line-of-sight direction (origin -> radar)
thr = theta_deg * pi / 180;
phr = phi_deg * pi / 180;
k_hat = [sin(thr) * cos(phr); sin(thr) * sin(phr); cos(thr)];

% Phase correction: exp(-j * 4*pi*f/c * r0.k_hat)
r0_dot_k = r0(1)*k_hat(1) + r0(2)*k_hat(2) + r0(3)*k_hat(3);
phase_correction = exp(-1j * (4 * pi / c) * r0_dot_k * freq_array);

Es_corrected = Es_f .* phase_correction;
end
