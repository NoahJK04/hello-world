function dxdt = rhs_psa_mol(~, x, par)
%RHS_PSA_MOL Right-hand side for PSA method-of-lines model.
%   State vector layout per bed: [P; yO2; Tg; Ts; wO2; wN2] (Nz each)
%   Then per bed voids: [nBtot;nBO2;TB;nTtot;nTO2;TT]
%   Receiver: [nRtot;nRO2;TR;TRs]

Nz = par.Nz;
R = par.R;

% Index helpers
idx = psa_indices(Nz);

% Unpack bed states
beds = cell(2,1);
for b = 1:2
    beds{b}.P = x(idx.bed(b).P);
    beds{b}.yO2 = x(idx.bed(b).yO2);
    beds{b}.Tg = x(idx.bed(b).Tg);
    beds{b}.Ts = x(idx.bed(b).Ts);
    beds{b}.wO2 = x(idx.bed(b).wO2);
    beds{b}.wN2 = x(idx.bed(b).wN2);
end

% Void and receiver states
voids = cell(2,1);
for b = 1:2
    voids{b}.nBtot = x(idx.void(b).nBtot);
    voids{b}.nBO2 = x(idx.void(b).nBO2);
    voids{b}.TB = x(idx.void(b).TB);
    voids{b}.nTtot = x(idx.void(b).nTtot);
    voids{b}.nTO2 = x(idx.void(b).nTO2);
    voids{b}.TT = x(idx.void(b).TT);
end

recv.nRtot = x(idx.recv.nRtot);
recv.nRO2 = x(idx.recv.nRO2);
recv.TR = x(idx.recv.TR);
recv.TRs = x(idx.recv.TRs);

% Pressures in voids/receiver
for b = 1:2
    voids{b}.PB = voids{b}.nBtot * R * voids{b}.TB / par.VB;
    voids{b}.PT = voids{b}.nTtot * R * voids{b}.TT / par.VT;
end
recv.PR = recv.nRtot * R * recv.TR / par.VR;

% Step information
step = par.step;

% Determine high/low bed indices
b_high = step.bed_high;
b_low = 3 - b_high;

% Flow computations between volumes (mol/s)
flows = struct();

% Feed to bottom void of high bed
val_feed = step.valves.feed;
flows.F_feed = valve_flow_linear(par.Pfeed, voids{b_high}.PB, val_feed, 0);

% Product from top void of high bed to receiver (flow controller)
val_prod = step.valves.product;
flows.F_prod = valve_flow_linear(voids{b_high}.PT, recv.PR, val_prod, step.flow_setpoints.product);

% Blowdown from bottom void of low bed to vent
val_blow = step.valves.blowdown;
flows.F_blow = valve_flow_linear(voids{b_low}.PB, par.Pamb, val_blow, 0);

% Purge from receiver to top void of low bed
val_purge = step.valves.purge;
flows.F_purge = valve_flow_linear(recv.PR, voids{b_low}.PT, val_purge, step.flow_setpoints.purge);

% Equalisation top and bottom (between beds)
val_eq_top = step.valves.eq_top;
val_eq_bottom = step.valves.eq_bottom;
flows.F_eq_top = valve_flow_linear(voids{b_high}.PT, voids{b_low}.PT, val_eq_top, 0);
flows.F_eq_bottom = valve_flow_linear(voids{b_high}.PB, voids{b_low}.PB, val_eq_bottom, 0);

% Bed through-flow estimate (positive bottom->top)
flow_bed = zeros(2,1);
flow_bed(b_high) = 0.5*(flows.F_feed + flows.F_prod);
flow_bed(b_low) = -0.5*(flows.F_purge + flows.F_blow);

% Initialize derivatives
for b = 1:2
    dP{b} = zeros(Nz,1);
    dyO2{b} = zeros(Nz,1);
    dTg{b} = zeros(Nz,1);
    dTs{b} = zeros(Nz,1);
    dwO2{b} = zeros(Nz,1);
    dwN2{b} = zeros(Nz,1);
end

dvoid = cell(2,1);
for b = 1:2
    dvoid{b}.nBtot = 0; dvoid{b}.nBO2 = 0; dvoid{b}.TB = 0;
    dvoid{b}.nTtot = 0; dvoid{b}.nTO2 = 0; dvoid{b}.TT = 0;
end

% Receiver derivatives

drecv.nRtot = 0; drecv.nRO2 = 0; drecv.TR = 0; drecv.TRs = 0;

% Bed PDE discretization
for b = 1:2
    P = beds{b}.P;
    yO2 = beds{b}.yO2;
    Tg = beds{b}.Tg;
    Ts = beds{b}.Ts;
    wO2 = beds{b}.wO2;
    wN2 = beds{b}.wN2;

    ctot = P ./ (R .* Tg);
    cO2 = yO2 .* ctot;
    yN2 = 1 - yO2;

    % Properties and dispersion
    EzO2 = zeros(Nz,1);
    EzN2 = zeros(Nz,1);
    rho_g = zeros(Nz,1);
    Cp_g = zeros(Nz,1);

    for i = 1:Nz
        props = gasprops(Tg(i), P(i), yO2(i), yN2(i));
        rho_g(i) = ctot(i) * (props.Mw/1000); % kg/m^3
        Cp_g(i) = props.Cp; % J/mol/K
        v = abs(flow_bed(b)) / max(ctot(i)*par.Ab, 1e-9);
        EzO2(i) = props.Dm_O2 + (0.73*v*par.rp / (1 + 9.49*props.Dm_O2/(max(v*par.rp,1e-9))));
        EzN2(i) = props.Dm_N2 + (0.73*v*par.rp / (1 + 9.49*props.Dm_N2/(max(v*par.rp,1e-9))));
    end

    % Inlet composition and temperature from connected void
    if flow_bed(b) >= 0
        cin_O2 = voids{b}.nBO2 / max(voids{b}.nBtot,1e-12) * (voids{b}.PB / (R*voids{b}.TB));
        Tin = voids{b}.TB;
    else
        cin_O2 = voids{b}.nTO2 / max(voids{b}.nTtot,1e-12) * (voids{b}.PT / (R*voids{b}.TT));
        Tin = voids{b}.TT;
    end

    % Convection terms using upwind
    v = flow_bed(b) ./ max(ctot .* par.Ab, 1e-9);
    v = mean(v); % use averaged velocity

    dcO2_dz = gradient_upwind(cO2, par.dz, v);
    dctot_dz = gradient_upwind(ctot, par.dz, v);

    % Dispersion term
    d2cO2_dz2 = second_derivative(cO2, par.dz, cin_O2);

    % Kinetics
    dwO2_vec = zeros(Nz,1);
    dwN2_vec = zeros(Nz,1);
    for i = 1:Nz
        [wO2_star, wN2_star] = iast_binary_sips(Tg(i), P(i)*par.u.Pa_to_bar, yO2(i), yN2(i), par.iso);
        dwO2_vec(i) = kinetics_qdf_darken(Tg(i), P(i)*par.u.Pa_to_bar, wO2(i), wO2_star, par.kin.O2, par.rp);
        dwN2_vec(i) = kinetics_qdf_darken(Tg(i), P(i)*par.u.Pa_to_bar, wN2(i), wN2_star, par.kin.N2, par.rp);
    end

    % Overall and component balances
    dctot_dt = -(1/par.eps_B) * (v .* dctot_dz) - (1-par.eps_B)/par.eps_B * par.rhos .* (dwO2_vec + dwN2_vec);
    dco2_dt = -(1/par.eps_B) * (v .* dcO2_dz) + (par.eps_i/par.eps_B) .* (EzO2 .* d2cO2_dz2) ...
        - (1-par.eps_B)/par.eps_B * par.rhos .* dwO2_vec;

    dy_dt = (dco2_dt - yO2 .* dctot_dt) ./ max(ctot,1e-9);

    % Gas temperature (simplified energy balance)
    hgs = 50; % W/m^2/K
    dTg_dt = -(v .* gradient_upwind(Tg, par.dz, v)) ...
        + (hgs*par.ap*(Ts - Tg)) ./ max(par.eps_B .* ctot .* Cp_g,1e-9) ...
        - (par.hw*(Tg - par.Tamb)) ./ max(par.eps_B .* ctot .* Cp_g,1e-9);

    % Solid temperature
    dTs_dt = (hgs*par.ap*(Tg - Ts)) ./ (par.rhos*par.Cps) ...
        + (-par.dH_O2.*dwO2_vec - par.dH_N2.*dwN2_vec) ./ (par.rhos*par.Cps);

    % Pressure derivative from ctot and Tg
    dP_dt = R .* (Tg .* dctot_dt + ctot .* dTg_dt);

    dP{b} = dP_dt;
    dyO2{b} = dy_dt;
    dTg{b} = dTg_dt;
    dTs{b} = dTs_dt;
    dwO2{b} = dwO2_vec;
    dwN2{b} = dwN2_vec;
end

% Void balances (ideal mixing)
for b = 1:2
    yB_O2 = voids{b}.nBO2 / max(voids{b}.nBtot,1e-12);
    yT_O2 = voids{b}.nTO2 / max(voids{b}.nTtot,1e-12);

    % Bottom void flows
    Fin_B = 0; Fout_B = 0; y_in_B = yB_O2; Tin_B = voids{b}.TB;
    if b == b_high
        Fin_B = flows.F_feed + flows.F_eq_bottom;
        y_in_B = par.yO2_feed;
        Tin_B = par.Tamb;
        Fout_B = max(flow_bed(b),0);
    else
        Fin_B = max(-flow_bed(b),0);
        Fout_B = flows.F_blow + flows.F_eq_bottom;
    end

    dvoid{b}.nBtot = Fin_B - Fout_B;
    dvoid{b}.nBO2 = Fin_B * y_in_B - Fout_B * yB_O2;
    dvoid{b}.TB = 0; % isothermal void approximation

    % Top void flows
    Fin_T = 0; Fout_T = 0; y_in_T = yT_O2; Tin_T = voids{b}.TT;
    if b == b_high
        Fin_T = max(flow_bed(b),0);
        Fout_T = flows.F_prod + flows.F_eq_top;
    else
        Fin_T = flows.F_purge + flows.F_eq_top;
        y_in_T = recv.nRO2 / max(recv.nRtot,1e-12);
        Tin_T = recv.TR;
        Fout_T = max(-flow_bed(b),0);
    end

    dvoid{b}.nTtot = Fin_T - Fout_T;
    dvoid{b}.nTO2 = Fin_T * y_in_T - Fout_T * yT_O2;
    dvoid{b}.TT = 0;
end

% Receiver balances
prod_y = voids{b_high}.nTO2 / max(voids{b_high}.nTtot,1e-12);

Fin_R = flows.F_prod;
Fout_R = flows.F_purge; % purge drawn from receiver

recv_y = recv.nRO2 / max(recv.nRtot,1e-12);

% Species
ndot_R = Fin_R - Fout_R;
ndot_RO2 = Fin_R * prod_y - Fout_R * recv_y;

% Energy (simple mixing + shell heat exchange)
props_R = gasprops(recv.TR, recv.PR, recv_y, 1-recv_y);
Cp_R = props_R.Cp;

TR_dot = 0;
if recv.nRtot > 1e-9
    TR_dot = (Fin_R*Cp_R*(par.Tamb - recv.TR) - Fout_R*Cp_R*(recv.TR - recv.TR)) / (recv.nRtot*Cp_R);
end

TRs_dot = (par.hR*par.AR*(recv.TR - recv.TRs) - par.hRa*par.AR*(recv.TRs - par.Tamb)) / (par.mR*par.CpR);

% Assign receiver derivatives

drecv.nRtot = ndot_R;
drecv.nRO2 = ndot_RO2;
drecv.TR = TR_dot;
drecv.TRs = TRs_dot;

% Pack derivatives

dxdt = zeros(size(x));
for b = 1:2
    dxdt(idx.bed(b).P) = dP{b};
    dxdt(idx.bed(b).yO2) = dyO2{b};
    dxdt(idx.bed(b).Tg) = dTg{b};
    dxdt(idx.bed(b).Ts) = dTs{b};
    dxdt(idx.bed(b).wO2) = dwO2{b};
    dxdt(idx.bed(b).wN2) = dwN2{b};

    dxdt(idx.void(b).nBtot) = dvoid{b}.nBtot;
    dxdt(idx.void(b).nBO2) = dvoid{b}.nBO2;
    dxdt(idx.void(b).TB) = dvoid{b}.TB;
    dxdt(idx.void(b).nTtot) = dvoid{b}.nTtot;
    dxdt(idx.void(b).nTO2) = dvoid{b}.nTO2;
    dxdt(idx.void(b).TT) = dvoid{b}.TT;
end

dxdt(idx.recv.nRtot) = drecv.nRtot;
dxdt(idx.recv.nRO2) = drecv.nRO2;
dxdt(idx.recv.TR) = drecv.TR;
dxdt(idx.recv.TRs) = drecv.TRs;

end

function idx = psa_indices(Nz)
% Build indices for state vector
base = 0;
for b = 1:2
    idx.bed(b).P = base + (1:Nz);
    base = base + Nz;
    idx.bed(b).yO2 = base + (1:Nz);
    base = base + Nz;
    idx.bed(b).Tg = base + (1:Nz);
    base = base + Nz;
    idx.bed(b).Ts = base + (1:Nz);
    base = base + Nz;
    idx.bed(b).wO2 = base + (1:Nz);
    base = base + Nz;
    idx.bed(b).wN2 = base + (1:Nz);
    base = base + Nz;

    idx.void(b).nBtot = base + 1; base = base + 1;
    idx.void(b).nBO2 = base + 1; base = base + 1;
    idx.void(b).TB = base + 1; base = base + 1;
    idx.void(b).nTtot = base + 1; base = base + 1;
    idx.void(b).nTO2 = base + 1; base = base + 1;
    idx.void(b).TT = base + 1; base = base + 1;
end

idx.recv.nRtot = base + 1; base = base + 1;
idx.recv.nRO2 = base + 1; base = base + 1;
idx.recv.TR = base + 1; base = base + 1;
idx.recv.TRs = base + 1; base = base + 1;
end

function g = gradient_upwind(f, dz, v)
% Upwind gradient based on velocity sign
n = numel(f);
g = zeros(n,1);
if v >= 0
    g(1) = (f(1) - f(1)) / dz;
    for i = 2:n
        g(i) = (f(i) - f(i-1)) / dz;
    end
else
    g(n) = (f(n) - f(n)) / dz;
    for i = 1:n-1
        g(i) = (f(i+1) - f(i)) / dz;
    end
end
end

function d2 = second_derivative(f, dz, fin)
% Second derivative with Danckwerts inlet, zero-gradient outlet.
n = numel(f);
d2 = zeros(n,1);
% Ghost point at inlet using specified inlet concentration
f0 = fin;
for i = 2:n-1
    d2(i) = (f(i+1) - 2*f(i) + f(i-1)) / dz^2;
end
% Inlet
    d2(1) = (f(2) - 2*f(1) + f0) / dz^2;
% Outlet zero gradient (Neumann)
    d2(n) = (f(n) - 2*f(n) + f(n-1)) / dz^2;
end
