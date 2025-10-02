function params = psa_default_parameters(overrides)
%PSA_DEFAULT_PARAMETERS Create default PSA model parameter structure.
%   PARAMS = PSA_DEFAULT_PARAMETERS() returns a structure containing the
%   default thermodynamic, kinetic, and geometric parameters used by the
%   PSA nitrogen purification model. Optional overrides can be supplied by
%   passing a structure whose fields replace the defaults.
%
%   Example
%   -------
%       params = psa_default_parameters(struct('feed_pressure', 7e5));
%
%   See also PSA_DEFAULT_SCHEDULE, PSA_SIMULATION.

params = struct();
params.R = 8.314;               % J/mol-K
params.T = 298;                 % K
params.bed_length = 1.0;        % m (not explicitly used but included for completeness)
params.bed_diameter = 0.15;     % m
params.porosity = 0.38;         % void fraction (-)
params.skeletal_density = 1200; % kg/m^3
params.mass_transfer_coeff = 0.12; % 1/s (LDF rate)
params.henry_const_N2 = 0.12;   % mol/kg-bar
params.henry_const_O2 = 0.45;   % mol/kg-bar
params.reference_pressure = 1.0e5; % Pa (1 bar)
params.feed_composition = [0.78, 0.22]; % mol fraction [N2, O2]
params.product_purge_composition = [0.99, 0.01];
params.feed_temperature = params.T;
params.feed_pressure = 6.5e5;   % Pa
params.low_pressure = 1.2e5;    % Pa
params.initial_pressure = 1.5e5;% Pa
params.column_cycles = 4;       % number of steps per cycle
params.void_volume = params.porosity * ...
    (pi*(params.bed_diameter/2)^2 * params.bed_length);
params.solid_mass = (1 - params.porosity) * ...
    (pi*(params.bed_diameter/2)^2 * params.bed_length) * params.skeletal_density;
params.initial_yN2 = params.feed_composition(1);
params.initial_yO2 = params.feed_composition(2);

if nargin > 0 && ~isempty(overrides)
    override_fields = fieldnames(overrides);
    for k = 1:numel(override_fields)
        params.(override_fields{k}) = overrides.(override_fields{k});
    end
end

params.void_volume = params.porosity * ...
    (pi*(params.bed_diameter/2)^2 * params.bed_length);
params.solid_mass = (1 - params.porosity) * ...
    (pi*(params.bed_diameter/2)^2 * params.bed_length) * params.skeletal_density;

end
