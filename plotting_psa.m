function plotting_psa(results, par)
%PLOTTING_PSA Generate key PSA plots.

u = par.u;
idx = psa_indices(par.Nz);

x = results.x;
t = results.t;

nRtot = x(:, idx.recv.nRtot);
nRO2 = x(:, idx.recv.nRO2);
TR = x(:, idx.recv.TR);

PR = nRtot .* par.R .* TR / par.VR;

figure('Name','Receiver Pressure');
plot(t, PR * u.Pa_to_bar, 'LineWidth',1.5);
xlabel('Time [s]'); ylabel('Receiver Pressure [bar]'); grid on;

figure('Name','Product O2 mole fraction');
plot(t, nRO2 ./ max(nRtot,1e-12), 'LineWidth',1.5);
xlabel('Time [s]'); ylabel('y_{O2}'); grid on;

% Axial profiles at final time for bed 1
P1 = x(end, idx.bed(1).P)';
yO2_1 = x(end, idx.bed(1).yO2)';
wO2_1 = x(end, idx.bed(1).wO2)';

figure('Name','Bed 1 profiles');
subplot(2,1,1);
plot(par.z, yO2_1, 'LineWidth',1.5);
xlabel('z [m]'); ylabel('y_{O2}'); grid on;
subplot(2,1,2);
plot(par.z, wO2_1, 'LineWidth',1.5);
xlabel('z [m]'); ylabel('w_{O2} [kmol/kg]'); grid on;

end

function idx = psa_indices(Nz)
base = 0;
for b = 1:2
    idx.bed(b).P = base + (1:Nz); base = base + Nz;
    idx.bed(b).yO2 = base + (1:Nz); base = base + Nz;
    idx.bed(b).Tg = base + (1:Nz); base = base + Nz;
    idx.bed(b).Ts = base + (1:Nz); base = base + Nz;
    idx.bed(b).wO2 = base + (1:Nz); base = base + Nz;
    idx.bed(b).wN2 = base + (1:Nz); base = base + Nz;

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
