function u = helpers_units()
%HELPERS_UNITS Unit conversion helpers and constants.
%   Returns a struct with common unit conversions used in the PSA model.

u.bar_to_Pa = 1e5;           % Pa per bar
u.Pa_to_bar = 1e-5;          % bar per Pa
u.L_to_m3 = 1e-3;            % m^3 per L
u.cm3_to_m3 = 1e-6;          % m^3 per cm^3
u.mm_to_m = 1e-3;            % m per mm
u.hr_to_s = 3600;            % s per hour

% Normal conditions for "m^3_n" conversions
u.Tn = 273.15;               % K
u.Pn = 1.01325e5;            % Pa

% Ideal gas constant
u.R = 8.314462618;           % J/mol/K

% Molar masses (kg/kmol)
u.Mw_O2 = 31.998;            % kg/kmol
u.Mw_N2 = 28.0134;           % kg/kmol

% Convert normal volumetric flow [m^3_n/h] to mol/s
u.m3nph_to_molps = @(Qn) (Qn / u.hr_to_s) * (u.Pn / (u.R * u.Tn));
% Convert mol/s to normal volumetric flow [m^3_n/h]
u.molps_to_m3nph = @(F) (F * u.R * u.Tn / u.Pn) * u.hr_to_s;

end
