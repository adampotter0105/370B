% ME370B: Modeling and Advanced Concepts
% Project 9 - Part 3
% Andy Huynh

% Plot electrical potential (V) vs. current density (kA/m^2)
% Plot power density (kW/m^2) vs. current density (kA/m^2)

% Make a second version of this plot showing the first-law efficiency (LHV basis)
% and power density vs. current density (just to emphasize the effect on efficiency).

T_range = 700+273.15:100:1200+273.15;
D_H2_H2O_range = [20.4354313003926 25.6763702034351 31.6113729970254 38.2614224613519 45.6463392373703 53.7849307201175]*1e-4;
D_O2_N2_range = [1.79781502714731 2.14875657320207 2.52769556793826 2.93419328331759 3.36785305546872 3.82831337780325]*1e-4;
i_cathode_range = [54.3448348129138 171.928233567501 446.954515296547 1000 1989.74465120697 3606.02653665849];
i_net_max_range = [18130 43755 47097 50400 53707 56913];

% Given
p         = 1e5; % Pa, pressure in GDL
F         = 96485; % C/kmol, Faraday's constant
ion_cond_range = [0.902649181138542 2.74284526614914 6.89597012853954 15 29.1360516184125 51.7159093192034]./((2*F)^2);
i_anode_range = 100*i_cathode_range; % A/m2, exchange current density 
L_YSZ     = 50*1e-6; % 50 um
L_GDL     = 5*1e-3; % 5 mm
R         = 8.3145; % J/kmol/K

for k = 1:length(T_range)
    T = T_range(k);
    D_H2_H2O = D_H2_H2O_range(k);
    D_O2_N2 = D_O2_N2_range(k);
    i_cathode = i_cathode_range(k);
    i_net_max = i_net_max_range(k);

    ion_cond = ion_cond_range(k);
    i_anode = i_anode_range(k);

    c = p/(R*T);
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
    mu_anode = chemPotentials(gas_anode)/1e3; % J/mol
    LHV_mass0 = LHV_mass(gas_anode);
    MW_anode = meanMolecularWeight(gas_anode);
    LHV_mole = LHV_mass0*MW_anode;

    set(gas_cathode, 'T', T, 'P', p, 'X', x_cathode);
    mu_cathode = chemPotentials(gas_cathode)/1e3;
    
    % 1st Pass (Eq)
    mu_anode_e_eq = 0; % reference
    
    mu_anode_H2_eq   = mu_anode(iH2);   % J/mol
    mu_anode_H2O_eq  = mu_anode(iH2O);
    mu_cathode_O2_eq = mu_cathode(iO2);
    
    mu_YSZ_anode_O_eq   = mu_anode_H2O_eq + 2*mu_anode_e_eq - mu_anode_H2_eq;
    mu_YSZ_cathode_O_eq = mu_YSZ_anode_O_eq;
    mu_cathode_e_eq     = 0.5*mu_YSZ_cathode_O_eq - 0.25*mu_cathode_O2_eq;
    Phi_eq              = 0.5/F*(mu_anode_H2_eq + 0.5*mu_cathode_O2_eq - mu_anode_H2O_eq);  % J/C, electrical potential
    
    % 2nd pass (net rate of reactiton)
    steps = 500;
    i_net = linspace(0, i_net_max, steps); % A/m2, specified current     %%% what is the maximum current density?  => when x_cathode_O2 becomes 0.  
    v     = i_net/(2*F); % mol/s/m2, reaction velocity
    
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
    mu_GDL_anode_H2O_diff = R*T*log(1+J_anode_H2O*L_GDL/(x_anode(iH2O)*c*D_H2_H2O)); 
    mu_anode_O = mu_YSZ_anode_O_eq + mu_GDL_anode_H2O_diff + ...
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
    power        = i_net .* Phi_term; % (W/m^2)
    if k == 6
        power_max    = max(power);
    end
    % Losses
    ohmic_loss = mu_YSZ_O_diff/F;
    cathode_loss = (0.5*mu_cathode_O2 + 2*mu_cathode_e - mu_cathode_O)/F;
    anode_loss =  (mu_anode_H2 + mu_anode_O - mu_anode_H2O - 2*mu_anode_e_term)/F;
    gdl_loss = (mu_GDL_anode_H2_diff + mu_GDL_anode_H2O_diff + 0.5*mu_GDL_cathode_O2_diff)/F;
    
    % mole fractions
    x_anode_H2   = x_anode(iH2)*(1-(J_anode_H2*L_GDL)/(x_anode(iH2)*c*D_H2_H2O));
    x_anode_H2O  = x_anode(iH2O)*(1+(J_anode_H2O*L_GDL)/(x_anode(iH2O)*c*D_H2_H2O));

    Phi_term_plot(k,:) = Phi_term;
    i_net_plot(k,:) = i_net;
    power_plot(k,:) = power;
    %LHV_eff_plot(k,:) = power/LHV_mole;
    LHV_eff_plot(k,:) = Phi_term/LHV_mole; % THIS NEEDS TO BE CORRECTED
end

figure(1)
plot(i_net_plot(1,:)/1e3, Phi_term_plot(1,:), 'b-')
hold on
plot(i_net_plot(2,:)/1e3, Phi_term_plot(2,:), 'b-')
plot(i_net_plot(3,:)/1e3, Phi_term_plot(3,:), 'b-')
plot(i_net_plot(4,:)/1e3, Phi_term_plot(4,:), 'b-')
plot(i_net_plot(5,:)/1e3, Phi_term_plot(5,:), 'b-')
plot(i_net_plot(6,:)/1e3, Phi_term_plot(6,:), 'b-')
plot(i_net_plot(1,:)/1e3, power_plot(1,:)/power_max, 'g-')
plot(i_net_plot(2,:)/1e3, power_plot(2,:)/power_max, 'g-')
plot(i_net_plot(3,:)/1e3, power_plot(3,:)/power_max, 'g-')
plot(i_net_plot(4,:)/1e3, power_plot(4,:)/power_max, 'g-')
plot(i_net_plot(5,:)/1e3, power_plot(5,:)/power_max, 'g-')
plot(i_net_plot(6,:)/1e3, power_plot(6,:)/power_max, 'g-')
text(i_net_plot(1,end)/1e3,power_plot(1,end)/power_max,'700ºC')
text(i_net_plot(2,end)/1e3,power_plot(2,end)/power_max,'800ºC')
text(i_net_plot(3,end)/1e3,power_plot(3,end)/power_max,'900ºC')
text(i_net_plot(4,end)/1e3,power_plot(4,end)/power_max,'1000ºC')
text(i_net_plot(5,end)/1e3,power_plot(5,end)/power_max,'1100ºC')
text(i_net_plot(6,end)/1e3,power_plot(6,end)/power_max,'1200ºC')
legend('Electrical Potential (V)','Normalized Power Density')
xlabel('Current Density (kA/m^2)')
hold off


figure(2)
plot(i_net_plot(1,:)/1e3, LHV_eff_plot(1,:), 'b-')
hold on
plot(i_net_plot(2,:)/1e3, LHV_eff_plot(2,:), 'b-')
plot(i_net_plot(3,:)/1e3, LHV_eff_plot(3,:), 'b-')
plot(i_net_plot(4,:)/1e3, LHV_eff_plot(4,:), 'b-')
plot(i_net_plot(5,:)/1e3, LHV_eff_plot(5,:), 'b-')
plot(i_net_plot(6,:)/1e3, LHV_eff_plot(6,:), 'b-')
% plot(i_net_plot(1,:)/1e3, power_plot(1,:)/power_max, 'g-')
% plot(i_net_plot(2,:)/1e3, power_plot(2,:)/power_max, 'g-')
% plot(i_net_plot(3,:)/1e3, power_plot(3,:)/power_max, 'g-')
% plot(i_net_plot(4,:)/1e3, power_plot(4,:)/power_max, 'g-')
% plot(i_net_plot(5,:)/1e3, power_plot(5,:)/power_max, 'g-')
% plot(i_net_plot(6,:)/1e3, power_plot(6,:)/power_max, 'g-')
% text(i_net_plot(1,end)/1e3,power_plot(1,end)/power_max,'700ºC')
% text(i_net_plot(2,end)/1e3,power_plot(2,end)/power_max,'800ºC')
% text(i_net_plot(3,end)/1e3,power_plot(3,end)/power_max,'900ºC')
% text(i_net_plot(4,end)/1e3,power_plot(4,end)/power_max,'1000ºC')
% text(i_net_plot(5,end)/1e3,power_plot(5,end)/power_max,'1100ºC')
% text(i_net_plot(6,end)/1e3,power_plot(6,end)/power_max,'1200ºC')
legend('LHV Efficiency','Normalized Power Density')
xlabel('Current Density (kA/m^2)')
hold off