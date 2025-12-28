function F = valve_flow_linear(P_up, P_down, valve, setpoint)
%VALVE_FLOW_LINEAR Linear valve/flow controller model.
%   P_up, P_down [Pa], valve struct with mode and Cv.
%   setpoint [mol/s] used for flow_controller mode.
%   Returns molar flow F [mol/s] (positive from up to down).

u = helpers_units();

if ~valve.open
    F = 0;
    return;
end

if ~isfield(valve, 'Cv')
    valve.Cv = 0;
end

switch valve.mode
    case 'ball'
        % Treat as very large Cv (minimal resistance)
        Cv = 1e3; % kmol/s/bar
        dP_bar = (P_up - P_down) * u.Pa_to_bar;
        F = Cv * dP_bar * 1e3; % kmol/s -> mol/s
    case 'control'
        Cv = valve.Cv; % kmol/s/bar
        dP_bar = (P_up - P_down) * u.Pa_to_bar;
        F = Cv * dP_bar * 1e3; % mol/s
    case 'flow_controller'
        % Set magnitude, enforce direction
        if P_up >= P_down
            F = setpoint;
        else
            F = 0;
        end
    otherwise
        error('Unknown valve mode: %s', valve.mode);
end
end
