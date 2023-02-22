% Set 6 Problem 3b

clc;

% Make a P vs. T diagram

% Provide access to support files via the Matlab path.
addpath 'Fundamental_Relation_Files' 
addpath 'Fundamental_Relation_Data'
addpath 'Setup_Files' 
addpath 'Property_Files' 

% Clean up and get ready to go.
clear all
format compact
fprintf('\n************************************************************\n')

% Set up the basic storage and load the FR files.
Setup_Props_i;

% Set which of the loaded species you want to work with.  You might want to
% change this around to see what the plots look like for other species.
ispecies = nH2;

Temps = linspace(Ttrip_i(ispecies), Tcrit_i(ispecies)-.3);
Pres = zeros(length(Temps),1);

for j = 1:length(Temps)
    [Sat_pressure, rho_liq_sat, rho_vap_sat] = Saturation_iT_NR(Temps(j));
    Pres(j,1) = Sat_pressure;
    disp(Temps(j));
end

figure(1)
plot(Temps,Pres/1e6,'LineWidth',2);
xlabel('Temperature (K)');
ylabel('Pressure (MPa)');
title('P-T Diagram for nH_2');
xlim([0 40]);

figure(2)
semilogy(-1./Temps,Pres/1e6,'LineWidth',2);
xlabel('Temperature (K)');
ylabel('Pressure (MPa)');
title('Duhring Diagram for nH_2');




