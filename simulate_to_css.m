function [results, history] = simulate_to_css(x0, par)
%SIMULATE_TO_CSS Run repeated cycles until cyclic steady state.

max_cycles = par.max_cycles;
css_tol = par.css_tol;

x_prev = x0;
results = [];

history.err = [];

for n = 1:max_cycles
    res = simulate_one_cycle(x_prev, par);
    results = res;

    x_end = res.x_end;
    err = max(abs((x_end - x_prev) ./ max(abs(x_prev), 1e-8)));
    history.err(end+1,1) = err; %#ok<AGROW>

    if err < css_tol
        break;
    end

    x_prev = x_end;
end

history.cycles = n;
end
