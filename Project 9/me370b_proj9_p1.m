
% Constants
R = 8.314; % Gas constant (J/mol*K)
F = 9.648533e4; % Faraday's Constant (C/mol)
kb = 1.38065e-23; % Boltzman Constant (J/K) 


% Button Cell Parameters
T = 1000 + 273; % Opeerating Temperature (K)
P_gas = 1e5; % Gas pressure (Pa)
K_ion = 15; %  Ionic Conductivity (S/m)
D_H2_water = 3.8378e-3; % Binary Diffusivity (m^2/s)
D_O2_N2 = 2.9417e-4; % Binary Diffusivity (m^2/s)
l_elec = 50e-6; % Width of electrolyte (m)
l_gdl = 5e-3; % Width of GDL
J_cath = 1000; % Exchange current density (A/m^2)
J_an = 100*J_cath; % Exchange current density (A/m^2)


% Cathode Gas (engineering air)
g_cath = GRI30;
set(g_cath, "T", T, "P", P_gas, "X", "N2:0.79, O2:0.21")

% Anode Gas (Humid H2)
g_an = GRI30;
set(g_cath, "T", T, "P", P_gas, "X", "H2:0.93, H2O:0.03")
