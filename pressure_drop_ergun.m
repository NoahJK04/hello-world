function dPdz = pressure_drop_ergun(v, rho, mu, par)
%PRESSURE_DROP_ERGUN Ergun equation for packed bed pressure gradient.
%   v [m/s] superficial velocity, rho [kg/m^3], mu [Pa*s]
%   Returns dP/dz [Pa/m].

rp = par.rp;
psi = par.psi;
eps = par.eps_i;

term1 = 150 * (1-eps)^2 * mu * v / (psi^2 * rp^2 * eps^3);
term2 = 1.75 * (1-eps) * rho * v * abs(v) / (psi * rp * eps^3);

dPdz = -(term1 + term2);
end
