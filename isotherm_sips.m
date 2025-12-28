function w = isotherm_sips(T, p_bar, iso)
%ISOTHERM_SIPS Pure-component Sips isotherm with temperature dependence.
%   T [K], p_bar [bar], iso struct with fields ws0, chi, b0, Q, n0, alpha.
%   Returns loading w [kmol/kg].

T0 = iso.T0;
nR = 8.314462618;
ws = iso.ws0 * exp(iso.chi * (1 - T0 ./ T));
b = iso.b0 * exp((iso.Q / nR) * (1./T - 1./T0)); % Q in J/mol
n = iso.n0 + iso.alpha * (T./T0 - 1);

arg = (b .* p_bar).^(1 ./ n);
w = ws .* (arg ./ (1 + arg));
end
