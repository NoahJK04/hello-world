function dPdz = ergun_drop(u, rho_g, params)
% Ergun pressure gradient (positive for flow in +z direction)

eps_b = params.eps_b;
mu = params.mu;
dp = params.dp;

A = 150 * mu * (1 - eps_b)^2 / (eps_b^3 * dp^2);
B = 1.75 * rho_g * (1 - eps_b) / (eps_b^3 * dp);

dPdz = A * u + B * u * abs(u);
end
