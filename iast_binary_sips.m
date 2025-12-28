function [wO2, wN2, wtot] = iast_binary_sips(T, P_bar, yO2, yN2, iso)
%IAST_BINARY_SIPS Binary IAST for Sips pure-component isotherms.
%   Returns wO2, wN2, wtot [kmol/kg].

if yO2 < 1e-12
    wO2 = 0; wN2 = isotherm_sips(T, P_bar, iso.N2); wtot = wN2; return;
elseif yN2 < 1e-12
    wN2 = 0; wO2 = isotherm_sips(T, P_bar, iso.O2); wtot = wO2; return;
end

% Root solve for fictitious pressures p0_O2 and p0_N2
p0_guess = [max(yO2*P_bar,1e-6); max(yN2*P_bar,1e-6)];

if exist('fsolve','file') == 2
    options = optimset('Display','off');
    sol = fsolve(@(x) iast_residuals(x, T, P_bar, yO2, yN2, iso), p0_guess, options);
else
    % Fallback to fminsearch if optimization toolbox is unavailable
    obj = @(x) sum(iast_residuals(x, T, P_bar, yO2, yN2, iso).^2);
    options = optimset('Display','off','MaxFunEvals',5000,'MaxIter',5000);
    sol = fminsearch(obj, p0_guess, options);
end

p0_O2 = max(sol(1), 1e-9);
p0_N2 = max(sol(2), 1e-9);

xO2 = yO2*P_bar / p0_O2;
xN2 = yN2*P_bar / p0_N2;

wO2_pure = isotherm_sips(T, p0_O2, iso.O2);
wN2_pure = isotherm_sips(T, p0_N2, iso.N2);

wtot = 1 / (xO2 / wO2_pure + xN2 / wN2_pure);
wO2 = xO2 * wtot;
wN2 = xN2 * wtot;

end

function F = iast_residuals(x, T, P_bar, yO2, yN2, iso)
% Equations for IAST:
% 1) spreading pressure equality
% 2) yO2*P/p0_O2 + yN2*P/p0_N2 = 1
p0_O2 = max(x(1), 1e-9);
p0_N2 = max(x(2), 1e-9);

pi_O2 = spreading_pressure(T, p0_O2, iso.O2);
pi_N2 = spreading_pressure(T, p0_N2, iso.N2);

F(1,1) = pi_O2 - pi_N2;
F(2,1) = yO2*P_bar/p0_O2 + yN2*P_bar/p0_N2 - 1;
end

function pi_val = spreading_pressure(T, p0, iso)
% Compute spreading pressure by quadrature: pi = integral_0^{p0} w(p)/p dp
integrand = @(p) isotherm_sips(T, p, iso) ./ max(p,1e-12);
pi_val = integral(integrand, 0, p0, 'RelTol',1e-6,'AbsTol',1e-10);
end
