function [i_guess, losses] = button_FC(V, T, p, xH2, xO2)
% Returns current flux and struct of losses of a differential fuel cell
% element for a given Votlage, temperature, pressure, input H2s mole
% fraction, and input O2 mole fraction.

%% Given Parameters
%T         = 1000+273.15; % K
%p         = 1e5; % Pa, pressure in GDL
F         = 96485; % C/kmol, Faraday's constant
%ion_cond  = 15/((2*F)^2); % S/m, ionic conductivity of YSZ
%i_anode   = 100*1000; % A/m2, exchange current density 
%i_cathode = 1000;
L_YSZ     = 50*1e-6; % 50 um
L_GDL     = 5*1e-3; % 5 mm
R         = 8.3145; % J/kmol/K
%D_H2_H2O  = 3.8378*1e-3; % m2/s
%D_O2_N2   = 2.9417*1e-4;
c         = p/(R*T);


% Dynamic Parameters
T_range = 700+273.15:100:1200+273.15;
D_H2_H2O_range = [20.4354313003926 25.6763702034351 31.6113729970254 38.2614224613519 45.6463392373703 53.7849307201175]*1e-4;
D_O2_N2_range = [1.79781502714731 2.14875657320207 2.52769556793826 2.93419328331759 3.36785305546872 3.82831337780325]*1e-4;
i_cathode_range = [54.3448348129138 171.928233567501 446.954515296547 1000 1989.74465120697 3606.02653665849];
i_net_max_range = [18130 43755 47097 50400 53707 56913];

ion_cond_range = [0.902649181138542 2.74284526614914 6.89597012853954 15 29.1360516184125 51.7159093192034]./((2*F)^2);
i_anode_range = 100*i_cathode_range; % A/m2, exchange current density 

T = T_range(k);
D_H2_H2O = interp1(T_range, D_H2_H2O_range, T);
D_O2_N2 = interp1(T_range, D_O2_N2_range, T);
i_cathode = interp1(T_range, i_cathode_range, T);
i_max = interp1(T_range, i_net_max_range, T);
ion_cond = interp1(T_range, ion_cond_range, T);
i_anode = interp1(T_range, i_anode_range, T);

gas  = Solution('GRI30.yaml');
iH2  = speciesIndex(gas, 'H2');
iH2O = speciesIndex(gas, 'H2O');
iO2  = speciesIndex(gas, 'O2');
iN2 = speciesIndex(gas, 'N2');
%MW = molecularWeights(gas); 

%% Gas input
nsp       = nSpecies(gas);
x_anode   = zeros(nsp, 1);
x_anode(iH2, 1)   = xH2; % mole fraction
x_anode(iH2O, 1)  = 1-xH2;
x_cathode = zeros(nsp, 1);
x_cathode(iO2, 1) = xO2; 
x_cathode(iN2, 1) = 1-xO2;

gas_anode  = Solution('GRI30.yaml');
gas_cathode  = Solution('GRI30.yaml');

set(gas_anode, 'T', T, 'P', p, 'X', x_anode);
mu_anode = chemPotentials(gas_anode)/1e3; % J/mol
set(gas_cathode, 'T', T, 'P', p, 'X', x_cathode);
mu_cathode = chemPotentials(gas_cathode)/1e3;

%% 1st Pass (Eq)
mu_anode_e_eq = 0; % reference

mu_anode_H2_eq   = mu_anode(iH2);   % J/mol
mu_anode_H2O_eq  = mu_anode(iH2O);
mu_cathode_O2_eq = mu_cathode(iO2);

mu_YSZ_anode_O_eq   = mu_anode_H2O_eq + 2*mu_anode_e_eq - mu_anode_H2_eq;
mu_YSZ_cathode_O_eq = mu_YSZ_anode_O_eq;
mu_cathode_e_eq     = 0.5*mu_YSZ_cathode_O_eq - 0.25*mu_cathode_O2_eq;
Phi_eq              = 0.5/F*(mu_anode_H2_eq + 0.5*mu_cathode_O2_eq - mu_anode_H2O_eq);  % J/C, electrical potential

% Check voltage request
if V < 0 || V > Phi_eq
    fprintf("button_FC Error: Voltage request out of bounds! \n")
end

% Loop NR around i_net
i_guess = i_max/2;

%NR Values
nr_tol = 1e-4;
max_it = 20;
it = 0;
di_guess = 1e-3; % how much to vary guess

% 2nd pass variables independent of i_guess
% Reaction rate at each node
R_anode_eq   = i_anode/(2*F); % A/m2, area-specific net reaction rate (reaction velocity)
R_cathode_eq = i_cathode/(2*F);

%% Start Newton Raphson Loop
while it < max_it

        it = it + 1;
        i_prev = i_guess;
        
        i_try = i_guess*[1-di_guess, 1, 1+di_guess]; % Put all three passees into an array
        
        % 2nd pass (net rate of reactiton)
        v = i_try/(2*F); % mol/s/m2, reaction velocity

        % flux
        J_anode_H2  = v;
        J_anode_H2O = v;
        %J_e   = 2*v;
        J_O   = v;
        J_cathode_O2  = 0.5*v;
        
        % EC potential of O ion in anode
        mu_anode_e_diff = 0; % very small 
        mu_GDL_anode_H2_diff  = -R*T*log(1-J_anode_H2*L_GDL/(x_anode(iH2)*c*D_H2_H2O)); % J/kmol
        mu_GDL_anode_H2O_diff = R*T*log(1+J_anode_H2O*L_GDL/(x_anode(iH2O)*c*D_H2_H2O)); 
        mu_anode_O = mu_YSZ_anode_O_eq + mu_GDL_anode_H2_diff + ...
            R*T*log(v/R_anode_eq+exp((mu_GDL_anode_H2O_diff+2*mu_anode_e_diff)/(R*T)));
        mu_anode_H2 = mu_anode_H2_eq - mu_GDL_anode_H2_diff;
        mu_anode_H2O = mu_anode_H2O_eq + mu_GDL_anode_H2O_diff;
        
        % EC potential of O ion across YSZ (find mu_cathode_o)
        mu_YSZ_O_diff = J_O*L_YSZ/ion_cond; 
        mu_cathode_O  = mu_anode_O + mu_YSZ_O_diff;
        
        % chemical potential of O2 in cathode
        x_cathode_O2  = 1-(1-x_cathode(iO2))*exp(J_cathode_O2*L_GDL/(c*D_O2_N2)); 
        mu_GDL_cathode_O2_diff = -R*T*log(x_cathode_O2/x_cathode(iO2));
        mu_cathode_O2 = mu_cathode_O2_eq - mu_GDL_cathode_O2_diff;
        
        % EC potential of electron in cathode
        mu_cathode_e = mu_cathode_e_eq + 0.25*mu_GDL_cathode_O2_diff + ...
            0.5*R*T*log(v/R_cathode_eq + exp((mu_cathode_O - mu_YSZ_cathode_O_eq)/(R*T)));
        mu_cathode_e_diff = 0; % very small
        mu_cathode_e_term = mu_cathode_e + mu_cathode_e_diff;
        
        % electrical potential between terminals
        mu_anode_e_term = mu_anode_e_eq; 
        Phi_term = -1/F*(mu_cathode_e_term - mu_anode_e_term);
        % power        = i_net .* Phi_term;
        % power_max    = max(power);
        
        % % Losses
        losses.ohmic_loss = mu_YSZ_O_diff/F;
        losses.cathode_loss = (0.5*mu_cathode_O2 + 2*mu_cathode_e - mu_cathode_O)/F;
        losses.anode_loss =  (mu_anode_H2 + mu_anode_O - mu_anode_H2O - 2*mu_anode_e_term)/F;
        losses.gdl_loss = (mu_GDL_anode_H2_diff + mu_GDL_anode_H2O_diff + 0.5*mu_GDL_cathode_O2_diff)/F;
        
        % % mole fractions
        % x_anode_H2   = x_anode(iH2)*(1-(J_anode_H2*L_GDL)/(x_anode(iH2)*c*D_H2_H2O));
        % x_anode_H2O  = x_anode(iH2O)*(1+(J_anode_H2O*L_GDL)/(x_anode(iH2O)*c*D_H2_H2O));
        %x_cathode_O2 = 1-(1-x_cathode(iO2))*exp(J_cathode_O2*L_GDL/(c*D_O2_N2)); 
        
        % Check for Convergence
        if abs(Phi_term(2)-V)/Phi_eq < nr_tol
            return
        end
        
        % Newton Raphson Step
        err_high = Phi_term(3) - V;
        err_low = Phi_term(1) - V;
        dedi = (err_high - err_low)/(2*di_guess*i_guess);
        i_err = Phi_term(2) - V;
        i_guess = i_guess - i_err/dedi;
        
        % Bisection if out of bounds
        if i_guess > i_max
            i_guess = (i_max-i_prev)/2 + i_prev;
        elseif i_guess < 0
            i_guess = i_prev/2;
        end

end

fprintf("button_FC Error: Failed to find current to matching current to requested voltage! \n")
