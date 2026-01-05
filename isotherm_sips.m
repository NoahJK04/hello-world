function [q, k] = isotherm_sips(P_A, P_B, params_A, params_B)
% Competitive Sips isotherm with simple pressure-dependent LDF rate

P_A = max(P_A, 1e-12);
P_B = max(P_B, 1e-12);

qsat = params_A.qsat;
bA = params_A.b;
nA = params_A.n;
bB = params_B.b;
nB = params_B.n;

termA = (bA .* P_A) .^ nA;
termB = (bB .* P_B) .^ nB;

q = qsat * (termA ./ (1 + termA + termB));

P_bar = (P_A + P_B) / 1e5;
k = params_A.k0 .* (1 + params_A.b_kin .* P_bar);
end
