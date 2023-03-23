function Phi_eq = button_phi_eq(T, p, xH2, xO2)
% Returns current flux and struct of losses of a differential fuel cell
% element for a given Votlage, temperature, pressure, input H2s mole
% fraction, and input O2 mole fraction.

% Given Parameters
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

a_np = 2.745*10^-4; % for nonpolar gas pairs
b_np = 1.823;
a_w = 3.64*10^-4; % for water with a nonpolar gas
b_w = 2.334;

K_0 = 15; % ionic conductivity at 1000C (S/m)
E_a = 96500; % activation energy (J/mol)
T_0 = 1000+273.15; % K

Ea_ecd = 100000; % activation energy for cathode exchange current density (J/mol)
ECD_cat_0 = 1000; % cathode exchange current density at 1000C (A/m2)
ECD_an_0 = ECD_cat_0*100;

% critical points
Tc_H2 = 33.18; % K
Pc_H2 = 13*100000/101325; % atm
M_H2 = .002016*1000; % g/mol
Tc_H2O = 373.95+273.15; % K
Pc_H2O = 217.7; % atm
M_H2O = .018015*1000; % g/mol
Tc_O2 = 154.58; % K
Pc_O2 = 50.0343*100000/101325; % atm
M_O2 = .0319988*1000; % g/mol
Tc_N2 = 126.2; % K
Pc_N2 = 34*100000/101325; % atm
M_N2 = .0280134*1000; % g/mol

P_atm = p/101325; % converts pressure to atm for diff calc

%D_H2_H2O = ((a_w * (T/sqrt(Tc_H2*Tc_H2O))^b_w) * (Pc_H2*Pc_H2O)^(1/3) * (Tc_H2*Tc_H2O)^(5/12) * (1/M_H2 + 1/M_H2O)^(1/2))/P_atm;
%D_O2_N2 = ((a_np * (T/sqrt(Tc_N2*Tc_O2))^b_np) * (Pc_N2*Pc_O2)^(1/3) * (Tc_N2*Tc_O2)^(5/12) * (1/M_N2 + 1/M_O2)^(1/2))/P_atm;
%i_cathode = ECD_cat_0*exp(-Ea_ecd/R*(1/T - 1/T_0));
% i_max = interp1(T_range, i_net_max_range, T);
%i_max = i_cathode*100;
%ion_cond = K_0*exp(-E_a/R*(1/T - 1/T_0));
%i_anode = ECD_an_0*exp(-Ea_ecd/R*(1/T - 1/T_0));

gas  = Solution('GRI30.yaml');
iH2  = speciesIndex(gas, 'H2');
iH2O = speciesIndex(gas, 'H2O');
iO2  = speciesIndex(gas, 'O2');
iN2 = speciesIndex(gas, 'N2');
%MW = molecularWeights(gas); 

% Gas input
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

% 1st Pass (Eq)
mu_anode_e_eq = 0; % reference

mu_anode_H2_eq   = mu_anode(iH2);   % J/mol
mu_anode_H2O_eq  = mu_anode(iH2O);
mu_cathode_O2_eq = mu_cathode(iO2);

mu_YSZ_anode_O_eq   = mu_anode_H2O_eq + 2*mu_anode_e_eq - mu_anode_H2_eq;
%mu_YSZ_cathode_O_eq = mu_YSZ_anode_O_eq;
%mu_cathode_e_eq     = 0.5*mu_YSZ_cathode_O_eq - 0.25*mu_cathode_O2_eq;
Phi_eq              = 0.5/F*(mu_anode_H2_eq + 0.5*mu_cathode_O2_eq - mu_anode_H2O_eq);  % J/C, electrical potential
