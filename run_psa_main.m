function run_psa_main(use_fast)
%RUN_PSA_MAIN Main entry point for twin-bed N2-PSA simulation.
%   Builds parameters, initializes state, runs CSS, prints KPIs, plots.

if nargin < 1
    use_fast = true;
end

par = params_reference_case();

if use_fast
    par = configure_fast_mode(par);
end

% Initial conditions
x0 = initial_state_psa(par);

% Run to CSS
[results, history] = simulate_to_css(x0, par);

% KPIs
kpi = kpis_from_results(results, par);
if isfield(par, 'fast_mode') && par.fast_mode
    fprintf('Fast mode enabled: KPIs are approximate and may be unstable.\n');
    kpi.purity = NaN;
    kpi.recovery = NaN;
    kpi.productivity = NaN;
    kpi.air_demand = NaN;
end

fprintf('CSS cycles: %d\n', history.cycles);
fprintf('Purity (N2): %.4f\n', kpi.purity);
fprintf('Recovery: %.4f\n', kpi.recovery);
fprintf('Productivity [m3n/h/m3]: %.4f\n', kpi.productivity);
fprintf('Air demand [m3n air / m3n N2]: %.4f\n', kpi.air_demand);

% Convergence plot
if par.enable_plots
    figure('Name','CSS Convergence');
    plot(1:numel(history.err), history.err, '-o');
    xlabel('Cycle'); ylabel('Max relative change'); grid on;
end

% Plots
if par.enable_plots
    plotting_psa(results, par);
end

end

function par = configure_fast_mode(par)
% Reduce model size and disable IAST for fast smoke testing.
par.fast_mode = true;
par.use_iast = false;
par.Nz = 5;
par.z = linspace(0, par.Hb, par.Nz)';
par.dz = par.z(2)-par.z(1);
par.max_cycles = 1;
par.css_tol = 1e-2;
par.rel_tol = 1e-3;
par.abs_tol = 1e-6;
par.enable_plots = false;

% Shorten step durations to speed up cycles
for i = 1:numel(par.cycle.steps)
    par.cycle.steps(i).duration = par.cycle.steps(i).duration * 0.05;
end
par.cycle.half_cycle_time = sum([par.cycle.steps.duration]);
end
