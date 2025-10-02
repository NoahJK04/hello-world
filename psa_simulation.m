function results = psa_simulation(varargin)
%PSA_SIMULATION Dynamic model of a single-column PSA unit for nitrogen purification.
%   RESULTS = PSA_SIMULATION() runs a dynamic simulation of a simplified
%   pressure swing adsorption (PSA) process that separates nitrogen from
%   air. The simulation models one adsorption column using a lumped
%   (well-mixed) approach with a linear driving force (LDF) adsorption
%   kinetic expression for nitrogen and oxygen. The process cycle includes
%   pressurization, adsorption, counter-current blowdown, and purge steps.
%
%   The function returns a RESULTS structure with the following fields:
%       time            - Simulation time vector (s)
%       states          - State history matrix [nN2, nO2, qN2, qO2, P]
%       step_indices    - Index of the active step for each time point
%       schedule        - Cycle schedule definition
%       params          - Model parameter structure
%       metrics         - Struct with cycle-averaged performance metrics
%
%   PSA_SIMULATION('plot', true) will generate diagnostic plots showing the
%   pressure profile, gas-phase compositions, and adsorbed loadings.
%
%   This model is intentionally simplified to emphasize the key dynamics of
%   PSA operation and to provide an approachable starting point for further
%   refinement. Assumptions include:
%       * The column is represented as a lumped, well-mixed volume.
%       * Adsorption equilibria follow a linear (Henry's law) isotherm.
%       * Kinetics follow a single LDF mass transfer coefficient.
%       * Heat effects are neglected (isothermal operation).
%       * Step pressures are enforced with a first-order approach.
%
%   Example
%   -------
%       results = psa_simulation('plot', true);
%       fprintf('Product purity: %.2f%%\n', 100*results.metrics.product_purity);
%       fprintf('Nitrogen recovery: %.2f%%\n', 100*results.metrics.n2_recovery);
%
%   See also ODE15S, TRAPZ

% Parse optional arguments
p = inputParser;
p.addParameter('plot', false, @(v)islogical(v) && isscalar(v));
p.addParameter('cycles', 8, @(v)isnumeric(v) && isscalar(v) && v>=1);
p.addParameter('params', struct(), @(v)isstruct(v));
p.addParameter('schedule', [], @(v)isstruct(v) || isempty(v));
p.parse(varargin{:});
opts = p.Results;

% Build default parameter set and schedule if not provided
params = psa_default_parameters(opts.params);
if isempty(opts.schedule)
    schedule = psa_default_schedule(params);
else
    schedule = opts.schedule;
end

% Run the requested number of cycles
[time, states, step_indices, metrics] = run_psa_cycles(params, schedule, opts.cycles);

results = struct('time', time, ...
                 'states', states, ...
                 'step_indices', step_indices, ...
                 'schedule', schedule, ...
                 'params', params, ...
                 'metrics', metrics);

if opts.plot
    plot_psa_diagnostics(results);
end

end

function [time, states, step_indices, metrics] = run_psa_cycles(params, schedule, n_cycles)
%RUN_PSA_CYCLES Integrate the PSA dynamics over multiple cycles.
if nargin < 3
    n_cycles = 6;
end

n_steps = numel(schedule);
time = [];
states = [];
step_indices = [];

% Initialise metrics
metrics = struct('feed_total', 0.0, ...
                 'feed_N2', 0.0, ...
                 'product_total', 0.0, ...
                 'product_N2', 0.0, ...
                 'waste_total', 0.0, ...
                 'waste_O2', 0.0, ...
                 'cycle_time', 0.0);

% Initial state at low pressure with feed composition
n_total0 = params.initial_pressure * params.void_volume / (params.R * params.T);
state0 = [params.initial_yN2 * n_total0; ...
          params.initial_yO2 * n_total0; ...
          params.henry_const_N2 * (params.initial_pressure/params.reference_pressure); ...
          params.henry_const_O2 * (params.initial_pressure/params.reference_pressure); ...
          params.initial_pressure];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
current_time = 0.0;

for cycle = 1:n_cycles
    for step_id = 1:n_steps
        step = schedule(step_id);
        tspan = [0, step.duration];
        odefun = @(t, x) psa_step_ode(t, x, params, step);
        [t_seg, x_seg] = ode15s(odefun, tspan, state0, options);

        % Evaluate stream information for metrics
        flows = evaluate_step_flows(x_seg, params, step);
        t_global = t_seg + current_time;

        % Accumulate time-series data (avoid duplicated timestamps)
        if isempty(time)
            time = t_global;
            states = x_seg;
            step_indices = repmat(step_id, size(t_seg));
        else
            time = [time; t_global(2:end)]; %#ok<AGROW>
            states = [states; x_seg(2:end, :)]; %#ok<AGROW>
            step_indices = [step_indices; repmat(step_id, numel(t_seg) - 1, 1)]; %#ok<AGROW>
        end

        % Integrate step metrics using trapezoidal rule
        if step.feed_stream
            metrics.feed_total = metrics.feed_total + trapz(t_global, flows.F_in);
            metrics.feed_N2 = metrics.feed_N2 + trapz(t_global, flows.F_in .* flows.y_in(:, 1));
        end

        if step.product_stream
            metrics.product_total = metrics.product_total + trapz(t_global, flows.F_out);
            metrics.product_N2 = metrics.product_N2 + trapz(t_global, flows.F_out .* flows.y_out(:, 1));
        end

        if step.waste_stream
            metrics.waste_total = metrics.waste_total + trapz(t_global, flows.F_out);
            metrics.waste_O2 = metrics.waste_O2 + trapz(t_global, flows.F_out .* flows.y_out(:, 2));
        end

        current_time = t_global(end);
        state0 = x_seg(end, :)';
    end
end

metrics.cycle_time = current_time / n_cycles;
metrics.product_purity = metrics.product_N2 / max(metrics.product_total, eps);
metrics.n2_recovery = metrics.product_N2 / max(metrics.feed_N2, eps);
metrics.waste_O2_fraction = metrics.waste_O2 / max(metrics.waste_total, eps);

end

function dxdt = psa_step_ode(~, x, params, step)
%PSA_STEP_ODE Column material balance for a single PSA step.
nN2 = x(1);
nO2 = x(2);
qN2 = x(3);
qO2 = x(4);
P = max(x(5), 1e3);

n_total = max(nN2 + nO2, 1e-12);
yN2 = nN2 / n_total;
yO2 = nO2 / n_total;

% Equilibrium loadings (Henry isotherm, using bar units for consistency)
bar_pressure = P / params.reference_pressure;
qN2_star = params.henry_const_N2 * yN2 * bar_pressure;
qO2_star = params.henry_const_O2 * yO2 * bar_pressure;

dqN2dt = params.mass_transfer_coeff * (qN2_star - qN2);
dqO2dt = params.mass_transfer_coeff * (qO2_star - qO2);

[Fin, Fout, yin] = compute_step_flows(step, params, P);

Fout = max(Fout, 0.0);
Fin = max(Fin, 0.0);

% Gas-phase material balances
solid_term_N2 = params.solid_mass * dqN2dt;
solid_term_O2 = params.solid_mass * dqO2dt;

dnN2dt = Fin * yin(1) - Fout * yN2 - solid_term_N2;
dnO2dt = Fin * yin(2) - Fout * yO2 - solid_term_O2;

% Pressure approach to target value
% Use first-order dynamics to reach set-point during the step
if isfield(step, 'P_target') && isfield(step, 'tau_P') && step.tau_P > 0
    dPdt = (step.P_target - P) / step.tau_P;
else
    dPdt = 0;
end

dxdt = [dnN2dt; dnO2dt; dqN2dt; dqO2dt; dPdt];

end

function flows = evaluate_step_flows(states, params, step)
%EVALUATE_STEP_FLOWS Compute inlet/outlet flow properties for diagnostics.
n_points = size(states, 1);
Fin = zeros(n_points, 1);
Fout = zeros(n_points, 1);
y_in = zeros(n_points, 2);
y_out = zeros(n_points, 2);

for i = 1:n_points
    nN2 = states(i, 1);
    nO2 = states(i, 2);
    P = states(i, 5);
    [Fin(i), Fout(i), yin] = compute_step_flows(step, params, P);
    y_in(i, :) = yin(:)';
    n_total = max(nN2 + nO2, 1e-12);
    y_out(i, :) = [nN2, nO2] / n_total;
end

flows = struct('F_in', Fin, ...
               'F_out', Fout, ...
               'y_in', y_in, ...
               'y_out', y_out);

end

function [Fin, Fout, yin] = compute_step_flows(step, params, pressure)
%COMPUTE_STEP_FLOWS Determine flowrates and compositions for each step.
yin = step.y_in;

switch step.type
    case 'pressurization'
        Fin = step.F_in;
        % Small leakage term to avoid pressure overshoot
        Fout = step.F_out + 0.02 * max(pressure - step.P_target, 0) / step.P_target * step.F_in;
    case 'adsorption'
        Fin = step.F_in;
        % Maintain near-equal inlet/outlet flows to keep pressure steady
        Fout = step.F_out + 0.01 * (pressure - step.P_target) / step.P_target * step.F_out;
    case 'depressurization'
        Fin = 0.0;
        deltaP = max(pressure - step.P_target, 0);
        Fout = step.F_out * (1 + 0.5 * deltaP / max(step.P_target, 1e3));
    case 'purge'
        Fin = step.F_in;
        Fout = step.F_out;
    otherwise
        Fin = step.F_in;
        Fout = step.F_out;
end

end

function plot_psa_diagnostics(results)
%PLOT_PSA_DIAGNOSTICS Generate diagnostic plots for the PSA simulation.
time = results.time / 60; % convert to minutes for readability
states = results.states;
step_indices = results.step_indices;
schedule = results.schedule;
metrics = results.metrics;

nN2 = states(:, 1);
nO2 = states(:, 2);
qN2 = states(:, 3);
qO2 = states(:, 4);
P = states(:, 5) / 1e5; % bar

yN2 = nN2 ./ max(nN2 + nO2, 1e-12);
yO2 = nO2 ./ max(nN2 + nO2, 1e-12);

figure('Name', 'PSA Nitrogen Purification Simulation', 'Color', 'w');

subplot(3, 1, 1);
plot(time, P, 'LineWidth', 1.5);
xlabel('Time (min)'); ylabel('Pressure (bar)');
title('Column Pressure Profile'); grid on;

subplot(3, 1, 2);
plot(time, yN2, 'LineWidth', 1.5); hold on;
plot(time, yO2, 'LineWidth', 1.2);
xlabel('Time (min)'); ylabel('Gas-phase mole fraction');
legend({'y_{N_2}', 'y_{O_2}'}, 'Location', 'best');
title('Gas-phase Composition'); grid on;

subplot(3, 1, 3);
plot(time, qN2, 'LineWidth', 1.5); hold on;
plot(time, qO2, 'LineWidth', 1.2);
xlabel('Time (min)'); ylabel('Loading (mol/kg)');
legend({'q_{N_2}', 'q_{O_2}'}, 'Location', 'best');
title('Adsorbed Phase Loading'); grid on;

sgtitle(sprintf(['PSA Cycle Diagnostics: Purity = %.1f%%, Recovery = %.1f%%', ...
                 '\nCycle time = %.1f min'], ...
                 100 * metrics.product_purity, ...
                 100 * metrics.n2_recovery, ...
                 metrics.cycle_time / 60));

% Annotate step transitions
hold on;
step_changes = [true; diff(step_indices) ~= 0];
transition_times = time(step_changes);
for k = 1:numel(transition_times)
    xline(transition_times(k), ':', 'Color', [0.6 0.6 0.6]);
end

% Add legend entries for steps
legend_entries = {schedule.name};
legend('Location', 'eastoutside'); %#ok<LLEG>
annotation('textbox', [0.78, 0.45, 0.2, 0.1], 'String', legend_entries, ...
           'FitBoxToText', 'on', 'EdgeColor', 'none');

end
