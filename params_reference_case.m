function par = params_reference_case()
%PARAMS_REFERENCE_CASE Build reference-case parameter struct for PSA model.
%   Returns struct par with geometry, isotherms, kinetics, valves, cycle,
%   gas properties, and numerical settings.

u = helpers_units();
par.u = u;

% -----------------------------
% Universal constants
% -----------------------------
par.R = u.R;                   % J/mol/K
par.Tamb = 293.15;             % K (20 C)
par.Pamb = 1.0 * u.bar_to_Pa;  % Pa
par.Pfeed = 8.0 * u.bar_to_Pa; % Pa (operating pressure)
par.yO2_feed = 0.209;
par.yN2_feed = 0.791;

% -----------------------------
% Geometry and bed properties
% -----------------------------
par.Db = 66 * u.mm_to_m;       % m
par.Hb = 0.581;                % m (average)
par.Ab = pi * (par.Db/2)^2;    % m^2
par.Vb = par.Ab * par.Hb;      % m^3

par.rhos = 711;                % kg/m^3
par.eps_i = 0.404;             % inter-particle voidage
par.eps_p = 0.234;             % intra-particle voidage
par.eps_B = par.eps_i + (1-par.eps_i)*par.eps_p;
par.rp = 0.830 * u.mm_to_m;    % m
par.ap = 3*(1-par.eps_i)/par.rp; % m^2/m^3 pellet (Eq. 3.1-7)
par.ksz = 0.675;               % W/m/K
par.Cps = 880;                 % J/kg/K
par.hw = 50;                   % W/m^2/K

% Void volumes
par.VB = 18.5 * u.cm3_to_m3;   % m^3 bottom void
par.VT = 35.0 * u.cm3_to_m3;   % m^3 top void

% Receiver
par.VR = 12 * u.L_to_m3;       % m^3
par.DR = 220 * u.mm_to_m;      % m
par.mR = 4.8;                  % kg
par.CpR = 500;                 % J/kg/K
par.hR = 50;                   % W/m^2/K
par.hRa = 20;                  % W/m^2/K
par.AR = pi * par.DR * (par.DR/2); % m^2 (approx. shell area)

% Adsorbent heats of adsorption (placeholder constants)
par.dH_O2 = -12e3;             % J/mol
par.dH_N2 = -10e3;             % J/mol

% -----------------------------
% Isotherm parameters (Sips)
% -----------------------------
par.iso.T0 = 293.15;
par.iso.O2.ws0 = 3.384e-3; % kmol/kg
par.iso.O2.chi = 1.104;
par.iso.O2.b0 = 9.436e-2;  % 1/bar
par.iso.O2.Q = 1.222e4;    % kJ/kmol (numerically equals J/mol)
par.iso.O2.n0 = 1.120;
par.iso.O2.alpha = 3.341e-1;
par.iso.O2.T0 = par.iso.T0;

par.iso.N2.ws0 = 2.707e-3; % kmol/kg
par.iso.N2.chi = 1.146;
par.iso.N2.b0 = 1.205e-1;  % 1/bar
par.iso.N2.Q = 1.187e4;    % kJ/kmol (numerically equals J/mol)
par.iso.N2.n0 = 1.185;
par.iso.N2.alpha = 2.263e-1;
par.iso.N2.T0 = par.iso.T0;

% -----------------------------
% Kinetics parameters
% -----------------------------
par.kin.O2.A = 1.208e-6;         % m^2/s
par.kin.O2.Ea = 1622.33;         % J/mol
par.kin.O2.A_over_rp2 = 1.7534;  % 1/s

par.kin.N2.A = 8.186e-5;         % m^2/s
par.kin.N2.Ea = 25098.55;        % J/mol
par.kin.N2.A_over_rp2 = 118.8268; % 1/s

% -----------------------------
% Momentum / Ergun parameters
% -----------------------------
par.psi = 1.959;               % shape factor

% -----------------------------
% Valve Cv parameters (kmol/s/bar)
% -----------------------------
par.valves.V1 = struct('mode','control','Cv',9.0e-4);
par.valves.V2 = struct('mode','control','Cv',5.0e-1);
par.valves.V3 = struct('mode','control','Cv',5.0e-1);
par.valves.V4 = struct('mode','ball','Cv',1.0); % product pipeline (ball)
par.valves.V5 = struct('mode','ball','Cv',1.0);
par.valves.V6 = struct('mode','ball','Cv',1.0); % tail gas
par.valves.V7 = struct('mode','ball','Cv',1.0);
par.valves.V8 = struct('mode','flow_controller','Cv',0.0);
par.valves.V9 = struct('mode','flow_controller','Cv',0.0);
par.valves.V10 = struct('mode','control','Cv',2.5e-5);
par.valves.V11 = struct('mode','control','Cv',2.5e-5);
par.valves.V12 = struct('mode','control','Cv',1.1e-4);
par.valves.V13 = struct('mode','control','Cv',1.1e-4);
par.valves.V14 = struct('mode','control','Cv',9.5e-4);
par.valves.V15 = struct('mode','control','Cv',9.5e-4);
par.valves.V16 = struct('mode','flow_controller','Cv',0.0);

% Purge and product setpoints
par.purge_factor = 0.40;
par.purge_flow_m3nph = 0.0485; % m^3_n/h
par.purge_flow = u.m3nph_to_molps(par.purge_flow_m3nph); % mol/s
par.product_flow = 0.05 * par.purge_flow; % placeholder, tuned by controller

% -----------------------------
% Numerical settings
% -----------------------------
par.Nz = 70;
par.z = linspace(0, par.Hb, par.Nz)';
par.dz = par.z(2)-par.z(1);

par.abs_tol = 1e-7;
par.rel_tol = 1e-5;
par.use_iast = true;
par.max_cycles = 50;
par.css_tol = 1e-4;
par.enable_plots = true;

% -----------------------------
% Cycle definition
% -----------------------------
par.cycle = psa_cycle_definition(par);

end
