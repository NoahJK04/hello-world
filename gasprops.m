function props = gasprops(T, P, yO2, yN2)
%GASPROPS Gas mixture properties: MW, Cp, viscosity, diffusion, conductivity.
%   T [K], P [Pa]. Returns struct props.

u = helpers_units();

% Molar masses
Mw_O2 = u.Mw_O2; % kg/kmol
Mw_N2 = u.Mw_N2;

% Mixture molar mass (kg/kmol)
props.Mw = yO2*Mw_O2 + yN2*Mw_N2;

% Cp (J/mol/K) - simple constant values
Cp_O2 = 29.4; % J/mol/K
Cp_N2 = 29.1;
props.Cp = yO2*Cp_O2 + yN2*Cp_N2;
props.Cv = props.Cp - u.R;

% Viscosity (Pa*s) - simple linear mixing
mu_O2 = 2.07e-5 * (T/300)^(0.7);
mu_N2 = 1.76e-5 * (T/300)^(0.7);
props.mu = yO2*mu_O2 + yN2*mu_N2;

% Diffusion coefficients (m^2/s) from Appendix forms (P in Pa)
P_bar = P * u.Pa_to_bar;
props.Dm_O2 = 4.81266e-8 * T^1.5 / (P_bar * T^(-0.291)) * 66.383537;
props.Dm_N2 = 4.81266e-8 * T^1.5 / (P_bar * T^(-0.219)) * 39.77204;

% Thermal conductivity (W/m/K)
props.kg = (7.793e-11*T + 2.456e-9)*yO2 + (6.914e-11*T + 4.118e-9)*yN2;

end
