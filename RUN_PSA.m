clc; clear; close all;

%% PSA MODEL: Isothermal, 1D, Axially Dispersed (Method of Lines)
% States: c_O2, c_N2, q_O2, q_N2
% Pressure from ideal gas: P = (c_O2 + c_N2) R T

%% ==================== PARAMETERS ====================
params = struct();

% Geometry
params.L = 6.5;               % m
params.D = 0.9;               % m
params.A = (pi*params.D^2)/4; % m^2

% Physical constants
params.R = 8.314;             % J/mol/K
params.T = 293;               % K

% Bed properties
params.eps_b = 0.404;          % bed void fraction
params.rho_s = 900;            % kg/m^3
params.dp = 2*0.00083;         % m

% Gas properties
params.mu = 1.81e-5;           % Pa*s
params.M_O2 = 0.032;           % kg/mol
params.M_N2 = 0.028;           % kg/mol

% Feed composition
params.y_feed_O2 = 0.21;
params.y_feed_N2 = 0.79;

% Pressures (Pa)
params.P_low = 1.0e5;          % Pa
params.P_high = 8.4e5;         % Pa
params.P_back = 7.5e5;         % Pa (adsorption backpressure)
params.P_atm = 1.01325e5;      % Pa

% Valve coefficients (mol/s/Pa^0.5)
params.Cv_in = 3.0e-2;
params.Cv_out = 2.0e-2;

% Axial grid
params.nz = 100;
params.z = linspace(0, params.L, params.nz)';
params.dz = params.z(2) - params.z(1);

% Sips isotherm parameters (competitive)
params.O2.qsat = 1.167;  % mol/kg
params.O2.b = 3.866e-3;  % 1/Pa
params.O2.n = 0.894;
params.O2.k0 = 5.3e-2;   % 1/s
params.O2.b_kin = 0.675; % 1/bar

params.N2.qsat = 0.681;
params.N2.b = 9.119e-3;
params.N2.n = 0.806;
params.N2.k0 = 2.4e-3;
params.N2.b_kin = 0.716;

%% ==================== TIME STEPS ====================
params.t_press = 6;   % s
params.t_ads = 20;    % s
params.t_blow = 5;    % s

%% ==================== INITIAL CONDITIONS ====================
P_init = params.P_low;
C_init = P_init/(params.R*params.T);

cO2_0 = 0.001*C_init*ones(params.nz,1);
cN2_0 = 0.999*C_init*ones(params.nz,1);

P_O2_init = cO2_0(1)*params.R*params.T;
P_N2_init = cN2_0(1)*params.R*params.T;

qO2_eq = isotherm_sips(P_O2_init, P_N2_init, params.O2, params.N2);
qN2_eq = isotherm_sips(P_N2_init, P_O2_init, params.N2, params.O2);

qO2_0 = 0.01*qO2_eq*ones(params.nz,1);
qN2_0 = 0.01*qN2_eq*ones(params.nz,1);

y0 = [cO2_0; cN2_0; qO2_0; qN2_0];

%% ==================== SOLVE PRESSURIZATION ====================
fprintf('Solving pressurization step...\n');
opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.05);
[t1, y1] = ode15s(@(t,y) rhs_psa(t,y,params,'press'), [0 params.t_press], y0, opts);

%% ==================== SOLVE ADSORPTION ====================
fprintf('Solving adsorption step...\n');
y_ads0 = y1(end,:)';
t_ads0 = params.t_press;
[t2, y2] = ode15s(@(t,y) rhs_psa(t,y,params,'ads'), [t_ads0 t_ads0+params.t_ads], y_ads0, opts);

%% ==================== SOLVE BLOWDOWN ====================
fprintf('Solving blowdown step...\n');
y_blow0 = y2(end,:)';
t_blow0 = t_ads0 + params.t_ads;
[t3, y3] = ode15s(@(t,y) rhs_psa(t,y,params,'blow'), [t_blow0 t_blow0+params.t_blow], y_blow0, opts);

%% ==================== POST PROCESS ====================
extract = @(Y) deal( ...
    Y(:,1:params.nz), ...
    Y(:,params.nz+1:2*params.nz), ...
    Y(:,2*params.nz+1:3*params.nz), ...
    Y(:,3*params.nz+1:4*params.nz));

[c1_O2, c1_N2, q1_O2, q1_N2] = extract(y1);
[c2_O2, c2_N2, q2_O2, q2_N2] = extract(y2);
[c3_O2, c3_N2, q3_O2, q3_N2] = extract(y3);

P1 = (c1_O2 + c1_N2)*params.R*params.T/1e5;
P2 = (c2_O2 + c2_N2)*params.R*params.T/1e5;
P3 = (c3_O2 + c3_N2)*params.R*params.T/1e5;

t_all = [t1; t2; t3];
cO2_all = [c1_O2; c2_O2; c3_O2];
cN2_all = [c1_N2; c2_N2; c3_N2];
qO2_all = [q1_O2; q2_O2; q3_O2];
qN2_all = [q1_N2; q2_N2; q3_N2];
P_all = [P1; P2; P3];

yO2_all = cO2_all ./ max(cO2_all + cN2_all, 1e-12);

%% ==================== DIAGNOSTICS ====================
fprintf('\nMass balance diagnostics (pressurization):\n');
mass_balance_check(t1, c1_O2, c1_N2, q1_O2, q1_N2, params, 'press');

fprintf('\nMass balance diagnostics (adsorption):\n');
mass_balance_check(t2, c2_O2, c2_N2, q2_O2, q2_N2, params, 'ads');

fprintf('\nMass balance diagnostics (blowdown):\n');
mass_balance_check(t3, c3_O2, c3_N2, q3_O2, q3_N2, params, 'blow');

%% ==================== PLOTS ====================
figure('Color','w');
subplot(1,2,1);
mesh(params.z, t_all, P_all);
xlabel('z (m)'); ylabel('t (s)'); zlabel('P (bar)');
title('Pressure'); view(-35,30); grid on;

subplot(1,2,2);
mesh(params.z, t_all, yO2_all);
xlabel('z (m)'); ylabel('t (s)'); zlabel('y_{O2}');
title('O_2 Mole Fraction'); view(-35,30); grid on;

figure('Color','w');
plot(t_all, P_all(:,end), 'LineWidth', 1.5);
hold on;
yline(params.P_high/1e5, 'r--');
yline(params.P_low/1e5, 'b--');
grid on;
xlabel('t (s)'); ylabel('P_{out} (bar)');
title('Outlet Pressure');

fprintf('\nSimulation complete.\n');
