function [results, history] = simulate_to_css(x0, par)
%SIMULATE_TO_CSS Run repeated cycles until cyclic steady state.

max_cycles = 50;
css_tol = 1e-4;

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
