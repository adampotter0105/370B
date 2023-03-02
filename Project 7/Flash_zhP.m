function [q, V, y, x, rg, rf, T] = Flash_zhP(z, h, P)
% Returns mass quality, vapor mole fraciton, vapor composition, liquid
% composition, vapor density, liquid density, and temperature for a given
% mixture of composition z under the vapor dome at enthalpy h and pressure P

% NR Parameters
h_tol = 0.001*h;
it = 0;
it_max = 20; % Max NR interations
Tinc = 0.05; % Temp increment for central finite differences

% [Thigh, ~, ~, ~] = Dew_cP(z,P);
% [Tlow, ~, ~, ~] = Bubble_cP(z,P);

% Initialize T and h
[~, isp] = max(z);
[r, T] = rT_ihP(isp, h, P);
h_guess = h_crT(z, r, T);

% Run NR to find T that matches h
while abs(h_guess-h)>h_tol
    Tlast = T;
    [q, V, y, x, rg, rf] = Flash_zTP(z,T,P); 
    % x-liq comp, y-vap comp, rg-vap density, rf-liq density, V-vapor mole
    % fraction, q-mass quality
    
    % Sum vap and liq enthalpies
    h_v = h_crT(y, rg, T);
    h_l = h_crT(x, rf, T);
    h_guess = h_v*V + (1-V)*h_l;

    % Central Finite Difference to find dh/dT for NR
    [~, V_low, y_low, x_low, rg_low, rf_low] = Flash_zTP(z,T-Tinc,P); 
    h_v_low = h_crT(y_low, rg_low, T-Tinc);
    h_l_low = h_crT(x_low, rf_low, T-Tinc);
    h_low = h_v_low*V_low + (1-V_low)*h_l_low;
    [~, V_high, y_high, x_high, rg_high, rf_high] = Flash_zTP(z,T+Tinc,P);
    h_v_high = h_crT(y_high, rg_high, T+Tinc);
    h_l_high = h_crT(x_high, rf_high, T+Tinc);
    h_high = h_v_high*V_high + (1-V_high)*h_l_high;
    dfdT = ((h_high-h)-(h_low-h))/(2*Tinc);
    
    % Use NR to calculate next Temp step
    dt = -(h_guess-h)/dfdT;
    T = T + dt;

    if T<Tlow
        T = 0.5*(Tlast-Tlow) + Tlow;
    elseif T>Thigh
        T = 0.5*(Thigh-Tlast) + Tlast;
    end
    
    if it<it_max
        fprintf("Flash_zhP failed to converge!! \n")
        break
    end

    it = it + 1;
end