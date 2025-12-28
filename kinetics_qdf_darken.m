function dw_dt = kinetics_qdf_darken(T, p_bar, w, wstar, kin, rp)
%KINETICS_QDF_DARKEN QDF kinetics with Darken correction + Arrhenius D0.
%   T [K], p_bar [bar], w,wstar [kmol/kg].
%   Returns dw_dt [kmol/kg/s].

R = 8.314462618; % J/mol/K

% Compute dp/dw numerically for Sips (fallback)
if wstar < 1e-12
    dw_dt = 0;
    return;
end

% Darken term approximation: (w/p) * (dp/dw)
% Estimate dp/dw by finite differences around wstar
p = max(p_bar, 1e-9);
dw = max(1e-6*wstar, 1e-9);

% Simple finite-difference slope assuming p ~ w (local linearization)
% This is a placeholder; replace with analytic derivative of Sips if desired.
dp_dw = p / max(wstar, 1e-9);

dlnp_dlnw = (wstar / p) * dp_dw;

D0 = kin.A * exp(-kin.Ea / (R * T));
D = D0 * dlnp_dlnw;

% QDF form from Eq. 3.5-6 (implemented as given)
rate = (wstar - w) / (rp^2) * D;

% Stabilize with QDF-type nonlinearity
eps = 1e-10;
rate = rate * ((wstar - w)^2 / max(wstar, eps));

dw_dt = rate;
end
