% Project 8 - Part 1
clear all
close all

% Given
T         = 1000+273.15; % K
p         = 101325; % Pa, pressure in GDL
ion_cond  = 15; % S/m, ionic conductivity of YSZ
F         = 96485*1e3; % C/kmol, Faraday's constant
i_anode   = 100*1000; % A/m2, exchange current density 
i_cathode = 1000;
L_YSZ     = 50*1e-6; % 50 um
L_GDL     = 5*1e-3; % 5 mm
R         = 8.3145*1e3; % J/kmol/K
D_H2_H2O  = 3.8378*1e-3; % m2/s
D_O2_N2   = 2.9417*1e-4;
c         = p/(R*T);

gas  = Solution('GRI30.yaml');
iH2  = speciesIndex(gas, 'H2');
iH2O = speciesIndex(gas, 'H2O');
iO2  = speciesIndex(gas, 'O2');
iN2 = speciesIndex(gas, 'N2');
MW = molecularWeights(gas); 

% gas input
nsp       = nSpecies(gas);
x_anode   = zeros(nsp, 1);
x_anode(iH2, 1)   = 0.97; % mole fraction
x_anode(iH2O, 1)  = 0.03;
x_cathode = zeros(nsp, 1);
x_cathode(iO2, 1) = 0.21; 
x_cathode(iN2, 1) = 0.79;

gas_anode  = Solution('GRI30.yaml');
gas_cathode  = Solution('GRI30.yaml');

set(gas_anode, 'T', T, 'P', p, 'X', x_anode);
mu_anode = chemPotentials(gas_anode);
set(gas_cathode, 'T', T, 'P', p, 'X', x_cathode);
mu_cathode = chemPotentials(gas_cathode);

% 1st Pass (Eq)
mu_anode_e_eq = 0; % reference

mu_anode_H2_eq   = mu_anode(iH2);   % J/kmol
mu_anode_H2O_eq  = mu_anode(iH2O);
mu_cathode_O2_eq = mu_cathode(iO2);

mu_YSZ_anode_O_eq   = mu_anode_H2O_eq + 2*mu_anode_e_eq - mu_anode_H2_eq;
mu_YSZ_cathode_O_eq = mu_YSZ_anode_O_eq;
mu_cathode_e_eq     = 0.5*mu_YSZ_cathode_O_eq - 0.25*mu_cathode_O2_eq;
Phi_eq              = 0.5/F*(mu_anode_H2_eq + 0.5*mu_cathode_O2_eq - mu_anode_H2O_eq);  % J/C, electrical potential

% 2nd pass (net rate of reactiton)
steps = 500;
i_net = linspace(0, 51232.7, steps); % A/m2, specified current     %%% what is the maximum current density?  => when x_cathode_O2 becomes 0.  
v     = i_net/(2*F); % kmol/s/m2, reaction velocity

% flux
J_anode_H2  = v;
J_anode_H2O = v;
J_e   = 2*v;
J_O   = v;
J_cathode_O2  = 0.5*v;

% Reaction rate at each node
R_anode_eq   = i_anode/(2*F); % A/m2, area-specific net reaction rate (reaction velocity)
R_cathode_eq = i_cathode/(2*F);

% EC potential of O ion in anode
mu_anode_e_diff = 0; % very small 
mu_GDL_anode_H2_diff  = -R*T*log(1-J_anode_H2*L_GDL/(x_anode(iH2)*c*D_H2_H2O)); % J/kmol
mu_GDL_anode_H2O_diff = R*T*log(1-J_anode_H2O*L_GDL/(x_anode(iH2O)*c*D_H2_H2O)); 
mu_anode_O = mu_YSZ_anode_O_eq + mu_GDL_anode_H2O_diff + R*T*log(v/R_anode_eq+exp((mu_GDL_anode_H2O_diff+2*mu_anode_e_diff)/(R*T)));

% EC potential of O ion in cathode
mu_YSZ_O_diff = J_O*L_YSZ/ion_cond; 
mu_cathode_O  = mu_anode_O + mu_YSZ_O_diff;

% C potential of O2 in cathode
x_cathode_O2  = 1-(1-x_cathode(iO2))*exp(J_cathode_O2*L_GDL/(c*D_O2_N2)); 
mu_GDL_cathode_O2_diff = -R*T*log(x_cathode_O2/x_cathode(iO2));
mu_cathode_O2 = mu_cathode_O2_eq - mu_GDL_cathode_O2_diff;

% EC potential of electron in cathode
mu_cathode_e = mu_cathode_e_eq + 0.25*mu_GDL_cathode_O2_diff + 0.5*R*T*log(v/R_cathode_eq + exp((mu_cathode_O - mu_YSZ_cathode_O_eq)/(R*T)));
mu_cathode_e_diff = 0; % very small
mu_cathode_e_term = mu_cathode_e + mu_cathode_e_diff;

% electrical potential between terminals
mu_anode_e_term = mu_anode_e_eq; 
Phi_term = -1/F*(mu_cathode_e_term - mu_anode_e_term);
power        = i_net .* Phi_term;
power_max    = max(power);

% mole fractions
x_anode_H2   = x_anode(iH2)*(1-(J_anode_H2*L_GDL)/(x_anode(iH2)*c*D_H2_H2O));
x_anode_H2O  = x_anode(iH2O)*(1+(J_anode_H2O*L_GDL)/(x_anode(iH2O)*c*D_H2_H2O));
x_cathode_O2 = 1-(1-x_cathode(iO2))*exp(J_cathode_O2*L_GDL/(c*D_O2_N2)); 


%% plot
figure()
plot(i_net/1e3, Phi_term, 'b-')
hold on
plot(i_net/1e3, x_anode_H2, 'r-')
plot(i_net/1e3, x_anode_H2O, 'm-')
plot(i_net/1e3, x_cathode_O2, 'k-')
plot(i_net/1e3, power/power_max, 'g-')
yline(0.97,'r--','LineWidth',2)
yline(0.21,'k--','LineWidth',2)
yline(0.03,'m--','LineWidth',2)
text(5,1.1,sprintf('Max. Power:  %.1f%% (kW)',power_max/1e3))
legend('Electric potential (V)','Power/Max Power','Anode Hydrogen','Anode Water','Cathode Oxygen')
xlabel('Current density [kA/m^2]')
title('YSZ: 1000 ^oC 1 bar')
xlim([0 60])
ylim([0 1.2])
plotfixer



