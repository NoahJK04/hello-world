function run_psa_main()
%RUN_PSA_MAIN Main entry point for twin-bed N2-PSA simulation.
%   Builds parameters, initializes state, runs CSS, prints KPIs, plots.

par = params_reference_case();

% Initial conditions
x0 = initial_state(par);

% Run to CSS
[results, history] = simulate_to_css(x0, par);

% KPIs
kpi = kpis_from_results(results, par);

fprintf('CSS cycles: %d\n', history.cycles);
fprintf('Purity (N2): %.4f\n', kpi.purity);
fprintf('Recovery: %.4f\n', kpi.recovery);
fprintf('Productivity [m3n/h/m3]: %.4f\n', kpi.productivity);
fprintf('Air demand [m3n air / m3n N2]: %.4f\n', kpi.air_demand);

% Convergence plot
figure('Name','CSS Convergence');
plot(1:numel(history.err), history.err, '-o');
xlabel('Cycle'); ylabel('Max relative change'); grid on;

% Plots
plotting_psa(results, par);

end

function x0 = initial_state(par)
% Build initial state vector
Nz = par.Nz;
idx = psa_indices(Nz);

n_state = idx.recv.TRs;
x0 = zeros(n_state,1);

% Uniform bed initial conditions
P0 = par.Pamb;
T0 = par.Tamb;
yO2_0 = par.yO2_feed;

% Equilibrium loading at initial conditions
[wO2_star, wN2_star] = iast_binary_sips(T0, P0*par.u.Pa_to_bar, yO2_0, 1-yO2_0, par.iso);

for b = 1:2
    x0(idx.bed(b).P) = P0;
    x0(idx.bed(b).yO2) = yO2_0;
    x0(idx.bed(b).Tg) = T0;
    x0(idx.bed(b).Ts) = T0;
    x0(idx.bed(b).wO2) = wO2_star;
    x0(idx.bed(b).wN2) = wN2_star;

    % Bottom void
    x0(idx.void(b).nBtot) = P0 * par.VB / (par.R*T0);
    x0(idx.void(b).nBO2) = yO2_0 * x0(idx.void(b).nBtot);
    x0(idx.void(b).TB) = T0;

    % Top void
    x0(idx.void(b).nTtot) = P0 * par.VT / (par.R*T0);
    x0(idx.void(b).nTO2) = yO2_0 * x0(idx.void(b).nTtot);
    x0(idx.void(b).TT) = T0;
end

% Receiver
x0(idx.recv.nRtot) = P0 * par.VR / (par.R*T0);
x0(idx.recv.nRO2) = yO2_0 * x0(idx.recv.nRtot);
x0(idx.recv.TR) = T0;
x0(idx.recv.TRs) = T0;

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
