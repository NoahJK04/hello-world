function dydt = rhs_psa(t, y, params, mode)
% RHS for PSA steps: press, ads, blow

nz = params.nz;
dz = params.dz;
R = params.R;
T = params.T;
eps_b = params.eps_b;
rho_s = params.rho_s;

cO2 = y(1:nz);
cN2 = y(nz+1:2*nz);
qO2 = y(2*nz+1:3*nz);
qN2 = y(3*nz+1:4*nz);

cO2 = max(cO2, 1e-12);
cN2 = max(cN2, 1e-12);
cT = cO2 + cN2;
P = cT * R * T;

P_O2 = cO2 * R * T;
P_N2 = cN2 * R * T;

% Kinetics
[qO2_star, kO2] = isotherm_sips(P_O2, P_N2, params.O2, params.N2);
[qN2_star, kN2] = isotherm_sips(P_N2, P_O2, params.N2, params.O2);

dqO2dt = kO2 .* (qO2_star - qO2);
dqN2dt = kN2 .* (qN2_star - qN2);

sink = rho_s * ((1 - eps_b) / eps_b) * (dqO2dt + dqN2dt);

% Boundary conditions for velocity
switch mode
    case 'press'
        Pin = params.P_high + (params.P_low - params.P_high) * exp(-0.8*t);
        Psource = Pin;
        Pback = P; %#ok<NASGU>
        inlet_open = true;
        outlet_open = false;
    case 'ads'
        Psource = params.P_high;
        Pback = params.P_back;
        inlet_open = true;
        outlet_open = true;
    case 'blow'
        Psource = params.P_low;
        Pback = params.P_atm;
        inlet_open = false;
        outlet_open = true;
    otherwise
        error('Unknown mode');
end

[u_face, u_node, P_face] = compute_uP(P, sink, params, inlet_open, outlet_open, Psource, Pback);

% Dispersion
DL_O2 = axial_dispersion(abs(u_node), P, params.O2, params);
DL_N2 = axial_dispersion(abs(u_node), P, params.N2, params);
DL_O2 = max(DL_O2, 1e-12);
DL_N2 = max(DL_N2, 1e-12);

% Feed concentration
C_feed = Psource/(R*T);
cf_O2 = params.y_feed_O2 * C_feed;
cf_N2 = params.y_feed_N2 * C_feed;

% Inlet Danckwerts
u_in = u_face(1);
if inlet_open
    ghost_O2 = cO2(1) - (dz*u_in/DL_O2(1))*(cO2(1) - cf_O2);
    ghost_N2 = cN2(1) - (dz*u_in/DL_N2(1))*(cN2(1) - cf_N2);
else
    ghost_O2 = cO2(1);
    ghost_N2 = cN2(1);
end

% Outlet: zero gradient
ghost_top_O2 = cO2(end-1);
ghost_top_N2 = cN2(end-1);

L_O2 = [ghost_O2; cO2(1:end-1)];
R_O2 = [cO2(2:end); ghost_top_O2];
L_N2 = [ghost_N2; cN2(1:end-1)];
R_N2 = [cN2(2:end); ghost_top_N2];

disp_O2 = DL_O2 .* (R_O2 - 2*cO2 + L_O2) / dz^2;
disp_N2 = DL_N2 .* (R_N2 - 2*cN2 + L_N2) / dz^2;

% Convective fluxes at faces (upwind)
J_O2 = zeros(nz+1,1);
J_N2 = zeros(nz+1,1);

% Inlet face
if inlet_open
    J_O2(1) = u_face(1) * cf_O2;
    J_N2(1) = u_face(1) * cf_N2;
else
    J_O2(1) = 0;
    J_N2(1) = 0;
end

for k = 2:nz
    if u_face(k) >= 0
        J_O2(k) = u_face(k) * cO2(k-1);
        J_N2(k) = u_face(k) * cN2(k-1);
    else
        J_O2(k) = u_face(k) * cO2(k);
        J_N2(k) = u_face(k) * cN2(k);
    end
end

% Outlet face
if outlet_open
    J_O2(end) = u_face(end) * cO2(end);
    J_N2(end) = u_face(end) * cN2(end);
else
    J_O2(end) = 0;
    J_N2(end) = 0;
end

dJdz_O2 = (J_O2(2:end) - J_O2(1:end-1)) / dz;
dJdz_N2 = (J_N2(2:end) - J_N2(1:end-1)) / dz;

accum_O2 = rho_s * ((1 - eps_b) / eps_b) * dqO2dt;
accum_N2 = rho_s * ((1 - eps_b) / eps_b) * dqN2dt;

dcdt_O2 = disp_O2 - dJdz_O2 - accum_O2;
dcdt_N2 = disp_N2 - dJdz_N2 - accum_N2;

dydt = [dcdt_O2; dcdt_N2; dqO2dt; dqN2dt];

end
