function mass_balance_check(t, cO2_hist, cN2_hist, qO2_hist, qN2_hist, params, mode)

z = params.z(:)';
dz = params.dz;
A = params.A;
eps_b = params.eps_b;
rho_s = params.rho_s;
R = params.R;
T = params.T;

nt = numel(t);
nz = params.nz;

Ng_O2 = eps_b*A*trapz(z, cO2_hist, 2);
Ng_N2 = eps_b*A*trapz(z, cN2_hist, 2);
Na_O2 = (1-eps_b)*rho_s*A*trapz(z, qO2_hist, 2);
Na_N2 = (1-eps_b)*rho_s*A*trapz(z, qN2_hist, 2);

N_O2 = Ng_O2 + Na_O2;
N_N2 = Ng_N2 + Na_N2;

dN_O2 = gradient(N_O2, t);
dN_N2 = gradient(N_N2, t);

Fin_O2 = zeros(nt,1);
Fin_N2 = zeros(nt,1);
Fout_O2 = zeros(nt,1);
Fout_N2 = zeros(nt,1);

for k = 1:nt
    cO2 = cO2_hist(k,:)';
    cN2 = cN2_hist(k,:)';
    qO2 = qO2_hist(k,:)';
    qN2 = qN2_hist(k,:)';
    cT = cO2 + cN2;
    P = cT * R * T;

    [qO2_star, kO2] = isotherm_sips(cO2*R*T, cN2*R*T, params.O2, params.N2);
    [qN2_star, kN2] = isotherm_sips(cN2*R*T, cO2*R*T, params.N2, params.O2);
    dqO2dt = kO2 .* (qO2_star - qO2);
    dqN2dt = kN2 .* (qN2_star - qN2);
    sink = rho_s * ((1-eps_b)/eps_b) * (dqO2dt + dqN2dt);

    switch mode
        case 'press'
            Psource = params.P_high + (params.P_low - params.P_high) * exp(-0.8*t(k));
            Pback = P(end);
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

    [u_face, ~, ~] = compute_uP(P, sink, params, inlet_open, outlet_open, Psource, Pback);
    DL_O2 = axial_dispersion(abs(u_face(1:end-1)), P, params.O2, params);
    DL_N2 = axial_dispersion(abs(u_face(1:end-1)), P, params.N2, params);

    dcO2dz = (cO2(2) - cO2(1))/dz;
    dcN2dz = (cN2(2) - cN2(1))/dz;

    Fin_O2(k) = A * (u_face(1)*cO2(1) - eps_b*DL_O2(1)*dcO2dz);
    Fin_N2(k) = A * (u_face(1)*cN2(1) - eps_b*DL_N2(1)*dcN2dz);

    dcO2dzL = (cO2(end) - cO2(end-1))/dz;
    dcN2dzL = (cN2(end) - cN2(end-1))/dz;

    Fout_O2(k) = A * (u_face(end)*cO2(end) - eps_b*DL_O2(end)*dcO2dzL);
    Fout_N2(k) = A * (u_face(end)*cN2(end) - eps_b*DL_N2(end)*dcN2dzL);
end

R_O2 = dN_O2 - (Fin_O2 - Fout_O2);
R_N2 = dN_N2 - (Fin_N2 - Fout_N2);
R_tot = (dN_O2 + dN_N2) - ((Fin_O2+Fin_N2) - (Fout_O2+Fout_N2));

fprintf('  O2 max |res|: %.3e mol/s\n', max(abs(R_O2)));
fprintf('  N2 max |res|: %.3e mol/s\n', max(abs(R_N2)));
fprintf('  Tot max |res|: %.3e mol/s\n', max(abs(R_tot)));

end
