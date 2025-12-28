function result = simulate_one_cycle(x0, par)
%SIMULATE_ONE_CYCLE Integrate model over one full PSA cycle.
%   Returns struct with time vector and state history.

cycle = par.cycle;

x = x0;
result.t = [];
result.x = [];

for half = 1:2
    for i = 1:numel(cycle.steps)
        step = cycle.steps(i);
        if step.duration <= 0
            continue;
        end
        step.bed_high = mod(step.bed_high + (half-1), 2) + 1;
        par.step = step;

        tspan = [0 step.duration];
        opts = odeset('RelTol', par.rel_tol, 'AbsTol', par.abs_tol);
        [t, xhist] = ode15s(@(t, x) rhs_psa_mol(t, x, par), tspan, x, opts);

        if isempty(result.t)
            result.t = t;
            result.x = xhist;
        else
            result.t = [result.t; result.t(end) + t];
            result.x = [result.x; xhist];
        end

        x = xhist(end,:)';
    end
end

result.x_end = x;
end
