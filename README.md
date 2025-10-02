# MATLAB PSA Nitrogen Purification Toolkit

This repository contains MATLAB tooling for exploring a single-column pressure swing adsorption (PSA) process that upgrades air to high-purity nitrogen. You can interact with the process in real time using an app-style visualisation or drive studies programmatically with the core simulation API.

## Interactive visual simulator

Launch the visual explorer to tune parameters live and watch the nitrogen purity respond:

```matlab
psa_visual_app;
```

Key features:

- **Real-time controls** – adjust feed/low pressure set-points, step durations, flow rates, and kinetic coefficients with spinners and switches.
- **Cycle visualisation** – the app shades each step in the pressure and composition plots so you can see where the product stream is drawn.
- **Performance dashboard** – purity, recovery, cycle time, and average product flow update automatically after every run.
- **Manual or auto update** – leave auto-update on for immediate feedback or switch to manual mode to stage multiple changes before running a new case.

Use the *Reset defaults* button to return to the reference operating point. The UI builds on the same underlying model as the scripted workflow, so you can move between approaches seamlessly.

## Scripted simulation

The function `psa_simulation` exposes the PSA model for scripted studies:

```matlab
results = psa_simulation('plot', true);
```

By default the routine simulates eight cycles and returns a structure with state trajectories, step indices, the active schedule, and key metrics.

### Optional arguments

```matlab
results = psa_simulation('cycles', 12, ...
                         'plot', true, ...
                         'params', struct('feed_pressure', 7.0e5), ...
                         'schedule', custom_schedule);
```

- **`cycles`** – number of PSA cycles to simulate (default: 8).
- **`plot`** – enable built-in diagnostic plots.
- **`params`** – structure of overrides applied to the default parameter set (pressure levels, kinetics, sorbent properties, etc.).
- **`schedule`** – array of step definitions in the format returned by `psa_default_schedule`.

### Helper builders

Two convenience functions make it easy to configure bespoke studies:

```matlab
params = psa_default_parameters(struct('feed_pressure', 7.5e5, ...
                                       'mass_transfer_coeff', 0.18));
schedule = psa_default_schedule(params);
```

You can mutate the returned structures before passing them to `psa_simulation` to experiment with custom steps, timings, and flow policies.

## File overview

- `psa_visual_app.m` – interactive app for tuning PSA parameters and visualising results in real time.
- `psa_simulation.m` – core dynamic PSA model with helper routines for integration and metrics.
- `psa_default_parameters.m` – utility that constructs the default parameter set (with optional overrides).
- `psa_default_schedule.m` – utility that builds the baseline PSA cycle definition.
- `index.html`, `styles.css`, `images/` – legacy files from the original GitHub Pages starter template.

## Extending the model

Ideas for future work include multi-column scheduling, richer mixture models (argon, CO₂), non-isothermal energy balances, and distributed bed models. Contributions and experimentation are welcome!
