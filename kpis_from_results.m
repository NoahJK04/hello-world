function kpi = kpis_from_results(results, par)
%KPIS_FROM_RESULTS Compute purity, recovery, productivity, and air demand.

u = par.u;

% Extract last cycle data
x = results.x;
t = results.t;

% Receiver outlet composition approximation
idx = psa_indices(par.Nz);

nRtot = x(:, idx.recv.nRtot);
nRO2 = x(:, idx.recv.nRO2);
yO2 = nRO2 ./ max(nRtot, 1e-12);
yN2 = 1 - yO2;

% Product purity (average over cycle)
kpi.purity = mean(yN2);

% Approximate product flow based on controller setpoint
Fprod = par.cycle.steps(4).flow_setpoints.product; % mol/s
kpi.product_flow_m3nph = u.molps_to_m3nph(Fprod);

% Feed flow based on feed valve and pressure drop (approx)
Ffeed = par.valves.V1.Cv * ((par.Pfeed - par.Pamb) * u.Pa_to_bar) * 1e3; % mol/s

% Recovery: N2 recovered / N2 in feed
kpi.recovery = (Fprod * kpi.purity) / (Ffeed * par.yN2_feed);

% Productivity: m3_n/h N2 per m3 adsorbent
kpi.productivity = kpi.product_flow_m3nph / par.Vb;

% Air demand (m3_n/h air per m3_n/h N2)
kpi.air_demand = (u.molps_to_m3nph(Ffeed)) / max(kpi.product_flow_m3nph, 1e-9);

% Sanitize outputs
fields = {'purity','recovery','productivity','air_demand','product_flow_m3nph'};
for i = 1:numel(fields)
    if ~isfinite(kpi.(fields{i}))
        kpi.(fields{i}) = NaN;
    end
end

end

function idx = psa_indices(Nz)
% Duplicate index helper to avoid coupling
base = 0;
for b = 1:2
    base = base + Nz; % P
    base = base + Nz; % yO2
    base = base + Nz; % Tg
    base = base + Nz; % Ts
    base = base + Nz; % wO2
    base = base + Nz; % wN2
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
