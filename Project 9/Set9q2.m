% Project 9 Part 2

clc;

% Constants
R = 8.314; % Gas constant (J/mol*K)
F = 9.648533e4; % Faraday's Constant (C/mol)
kb = 1.38065e-23; % Boltzman Constant (J/K) 
h = 6.626e-34; % Plank's Constant

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
M_H2 = .002016; % kg/mol
Tc_H2O = 373.95+273.15; % K
Pc_H2O = 217.7; % atm
M_H2O = .018015; % kg/mol
Tc_O2 = 154.58; % K
Pc_O2 = 50.0343*100000/101325; % atm
M_O2 = .0319988; % kg/mol
Tc_N2 = 126.2; % K
Pc_N2 = 34*100000/101325; % atm
M_N2 = .0280134; % kg/mol

% temperature range
T = linspace(600+273.15,1300+273.15); % K

% results arrays
D_N2O2 = zeros(1, length(T)); 
D_H2H2O = zeros(1,length(T));
K = zeros(1, length(T));
ECD_cat = zeros(1,length(T));
ECD_an = zeros(1,length(T));

P = 100000; % Pa
P_atm = P/101325; % atm

cathode = GRI30();
set(cathode,'P',P, 'X', 'N2:0.79, O2:0.21');

anode = GRI30();
set(anode, 'P', P, 'X', 'H2:0.93, H2O:0.03');

for i = 1:length(T)
    set(cathode,'T',T(i));
    set(anode,'T',T(i));
    
    % Bird Stewart Lightfoot equation (P: atm, T: K, D: cm2/s)
    % cathode (N2, O2)
    D_N2O2(i) = ((a_np * (T(i)/sqrt(Tc_N2*Tc_O2))^b_np) * (Pc_N2*Pc_O2)^(1/3) * (Tc_N2*Tc_O2)^(5/12) * (1/M_N2 + 1/M_O2)^(1/2))/P_atm;
    % anode (H2, H2O)
    D_H2H2O(i) = ((a_w * (T(i)/sqrt(Tc_H2*Tc_H2O))^b_w) * (Pc_H2*Pc_H2O)^(1/3) * (Tc_H2*Tc_H2O)^(5/12) * (1/M_H2 + 1/M_H2O)^(1/2))/P_atm;
    
    % ionic conductivity
    K(i) = K_0*exp(-E_a/R*(1/T(i) - 1/T_0));
    
    % exchange current density
    ECD_cat(i) = ECD_cat_0*exp(-Ea_ecd/R*(1/T(i) - 1/T_0));
    ECD_an(i) = ECD_an_0*exp(-Ea_ecd/R*(1/T(i) - 1/T_0));
    
    
end

% get values to normalize
D_N2O2_1000 = ((a_np * (1273/sqrt(Tc_N2*Tc_O2))^b_np) * (Pc_N2*Pc_O2)^(1/3) * (Tc_N2*Tc_O2)^(5/12) * (1/M_N2 + 1/M_O2)^(1/2))/P_atm;
D_H2H2O_1000 = ((a_w * (1273/sqrt(Tc_H2*Tc_H2O))^b_w) * (Pc_H2*Pc_H2O)^(1/3) * (Tc_H2*Tc_H2O)^(5/12) * (1/M_H2 + 1/M_H2O)^(1/2))/P_atm;

% normalize
D_N2O2 = D_N2O2./D_N2O2_1000;
D_H2H2O = D_H2H2O./D_H2H2O_1000;
K = K./K_0;
ECD_cat = ECD_cat./ECD_cat_0;
ECD_an = ECD_an./ECD_an_0;

% plotting
T_plot = T - 273.15;

figure(1)
plot(T_plot,D_N2O2);
hold on;
plot(T_plot,D_H2H2O);
plot(T_plot,K);
plot(T_plot,ECD_cat);
% plot(T_plot,ECD_an);
legend('D_{N_2O_2}','D_{H_2H_2O}','Ionic Conductivity','Cathode Exchange Current Density');
xlabel(['Temperature (' char(176) 'C)']);
ylabel('Normalized Quantities');
hold off;





    
    