function DL = axial_dispersion(u, P, params_A, params)
% Simple correlation for axial dispersion

rp = params.dp/2;
eps_i = 0.404;

u = max(u, 1e-12);
P = max(P, 1e-12) / 1e5;

M1 = 28;
M2 = 32;
C = 3.5;
E = -0.25;
T = params.T;

omega = C * (T^E);
sigma = 3.6;

term1 = 1.86e-3 * (T^(3/2));
term2 = sqrt((1/M1) + (1/M2));
P_atm = P/1.01325;
term3 = P_atm .* (sigma^2) .* omega;
Dm = ((term1*term2) ./ term3) / 10000;

DL = 0.73.*Dm + (u.*rp) ./ ( eps_i .* (1 + 9.49.*((eps_i.*Dm) ./ (2.*u.*rp))) );
end
