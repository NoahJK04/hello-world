function Mmix = mix_mw(yO2, yN2, params)
% Mixture molecular weight

Mmix = yO2 * params.M_O2 + yN2 * params.M_N2;
end
