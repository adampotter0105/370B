clear all

% Constants
R = 8.314; % Gas constant (J/mol*K)
F = 9.648533e7; % Faraday's Constant (C/mol)
kb = 1.38065e-23; % Boltzman Constant (J/K) 
h = 6.626e-34; % Plank's Constant

% Button Cell Parameters
T = 1000 + 273; % Operating Temperature (K)
P = 1e5; % Gas pressure (Pa)
K_ion = 15; %  Ionic Conductivity (S/m)
D_H2_water = 3.8378e-3; % Binary Diffusivity (m^2/s)
D_O2_N2 = 2.9417e-4; % Binary Diffusivity (m^2/s)
l_elec = 50e-6; % Width of electrolyte (m)
l_gdl = 5e-3; % Width of GDL
J_cath = 1000; % Exchange current density (A/m^2)
J_anod = 100*J_cath; % Exchange current density (A/m^2)

% Cathode Gas (engineering air)
g_cath = GRI30;
set(g_cath, "T", T, "P", P, "X", 'N2:0.79,O2:0.21')

% Anode Gas (Humid H2)
g_anod = GRI30;
set(g_anod, "T", T, "P", P, "X", 'H2:0.97,H2O:0.03')

iO2 = speciesIndex(g_cath, 'O2');
iN2 = speciesIndex(g_cath, 'N2');
iH2 = speciesIndex(g_anod, 'H2');
iH2O = speciesIndex(g_anod, 'H2O');

% Electrochemical Potential    mu_elec = mu + x*F*(phi-phi_ref)

%% First Pass
mu_c = chemPotentials(g_cath);
mu_a = chemPotentials(g_anod);
phi_ocp = (1/(2*F))*(mu_a(iH2) + 0.5*mu_c(iO2) - mu_a(iH2O) )

% NEED mu_e_eq_a
% mu_e_eq_a = mu_e_eq_c + F*phi_ocp    (this requires mu_e_eq_c though)
mu_o_ysz_eq_a = mu_a(iH2O) + 2*mu_e_eq_a - mu_a(iH2);
mu_o_ysz_eq_c = mu_o_ysz_eq_a;

% Electron at Cathode
mu_e_eq_c = 0.5*mu_o_ysz_eq_c - 0.25*mu_c(iO2);

%% Second Pass
I = 1e3; % (A/m^2)  ASSIGNED CURRENT DENSITY
v = I/(2*F); % (mol/s*m^2)

% Diffusion Losses of GDL
c = P/(R*T);
dmu_gdl_h2 = -R*T*log(1-(J_anod*l_gdl)/(moleFraction(g_anod,'H2')*c*D_H2_water));
dmu_gdl_h2o = -R*T*log(1+(J_anod*l_gdl)/(moleFraction(g_anod,'H2O')*c*D_H2_water));
dmu_gdl_o2 = -R*T*log(( 1-(1-moleFraction(g_cath,'O2'))*exp(J_cath*l_gdl/ ...
    (c*D_O2_N2)) )/moleFraction(g_cath,'O2'));

% Adjust chem potentials based on GDL losses
mu_h2_a = mu_a(iH2) - dmu_gdl_h2;
mu_h2o_a = mu_a(iH2O) + dmu_gdl_h2o;
mu_o2_a = mu_a(iO2) - dmu_gdl_o2;

% Anode
% NEEDS R_eq_a
dmu_e_a = 0; % because Ni losses are zero, CHECK!
mu_o_a = mu_o_ysz_a + dmu_gdl_h2 + R*T*log(v/R_eq_a + ...
    exp((dmu_gdl_h2o+2*dmu_e_a)/(R*T)) );

% Ohmic Losses in YSZ Electrolyte 
dmu_ysz_o = J*l_elec/K_ion;
mu_o_c = mu_o_a +dmu_ysz_o;

% Cathode
% NEED R_eq_c
mu_e_c = mu_e_eq_c + 0.25*dmu_gdl_o2 + 0.5*R*T*ln(v/R_eq_c + exp((mu_o_c - ...
    mu_o_ysz_eq_c)/(R*T)) );

% Potential Between Electrodes (nickel losses neglected)
% NEED mu_e_a
dphi_ni = -(1/F)*(mu_e_c - mu_e_a);

