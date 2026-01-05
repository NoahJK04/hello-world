function [u_face, u_node, P_face] = compute_uP(P, sink, params, inlet_open, outlet_open, Psource, Pback)
% Compute u(z) from continuity and valve BCs, then P via Ergun.

nz = params.nz;
dz = params.dz;
eps_b = params.eps_b;
R = params.R;
T = params.T;

C = max(P, 1e-12)/(R*T);

% Inlet/outlet valve flows
if inlet_open
    Fin = params.Cv_in * sqrt(max(Psource - P(1), 0));
    Cin = Psource/(R*T);
    u_in = Fin / (params.A * max(Cin, 1e-12));
else
    u_in = 0;
end

if outlet_open
    Fout = params.Cv_out * sqrt(max(P(end) - Pback, 0));
    u_out = Fout / (params.A * max(C(end), 1e-12));
else
    u_out = 0;
end

% Integrate continuity: d(uC)/dz = -sink
flux = zeros(nz,1);
flux(1) = u_in * C(1);
for i = 2:nz
    loss = 0.5*(sink(i) + sink(i-1)) * dz;
    flux(i) = flux(i-1) - loss;
end

u_node = flux ./ max(C, 1e-12);

% Enforce outlet flow by scaling if needed
if outlet_open
    if u_node(end) > 1e-12
        scale = u_out / u_node(end);
        u_node = u_node * scale;
        flux = flux * scale;
    end
else
    u_node(end) = 0;
end

% Face velocities (average)
u_face = zeros(nz+1,1);
u_face(1) = u_in;
u_face(end) = u_out;
for i = 2:nz
    u_face(i) = 0.5*(u_node(i-1) + u_node(i));
end

% Pressure drop from Ergun using face velocities
P_face = zeros(nz+1,1);
P_face(1) = Psource;
for i = 2:nz+1
    u_f = u_face(i);
    P_guess = P_face(i-1);
    rho_g = max(P_guess, 1e-6) * mix_mw(params.y_feed_O2, params.y_feed_N2, params) / (R*T);
    dPdz = ergun_drop(u_f, rho_g, params);
    P_face(i) = max(P_face(i-1) - dPdz*dz, 1e-6);
end

end
