% Simple-cycle gas turbine burning natural gas.
% C.F. Edwards, 12/19/10
% Dongwon Ka 01/28/2023

clear all
format compact
fprintf('\n****************************************************************\n')

% Set the number of "differential" steps for the polytropic process.
steps = 200;

% Cycle specifications:
global To;
To = 20+273.15;
global Po;
Po = oneatm;
global mu_o

Pressure_Ratio = 23;   % pressure ratio 23:1
Burner_Pressure_Ratio                 = 0.95;   %%%
Air_Compressor_Polytropic_Efficiency  = 0.90;
Fuel_Compressor_Polytropic_Efficiency = 0.90;
Mix_Turbine_Polytropic_Efficiency     = 0.88;
Tmax = 1700;

% Make a gas object to work with in Cantera.
gas   = Solution('gri30.yaml');
N     = nSpecies(gas);
iCH4  = speciesIndex(gas,'CH4');
iC2H6 = speciesIndex(gas,'C2H6');
iC3H8 = speciesIndex(gas,'C3H8');
iCO2  = speciesIndex(gas,'CO2');
iO2   = speciesIndex(gas,'O2');
iN2   = speciesIndex(gas,'N2');
iAR   = speciesIndex(gas, 'AR');
iH2O  = speciesIndex(gas, 'H2O');
M     = molecularWeights(gas);


% Start with RH=50% air at ambient conditions.
mdot_air = 1;
water=Solution('liquidvapor.yaml','water');
P_sat_water=satPressure(water,To);
RH=0.5; % 50%
P_H2O=RH*P_sat_water;
xair = zeros(1,N);
xair(iH2O)=P_H2O/Po;
xair(iN2) = (1-xair(iH2O))*0.757223/(1-0.031208);
xair(iO2) = (1-xair(iH2O))*0.202157/(1-0.031208);
xair(iAR) = (1-xair(iH2O))*0.009015/(1-0.031208);
xair(iCO2)= (1-xair(iH2O))*0.000397/(1-0.031208);

Ta1 = To;
Pa1 = Po;
set(gas,'T',Ta1,'P',Pa1,'X',xair);
mu_o = chemPotentials(gas);
M_air = meanMolecularWeight(gas);
sa1 = entropy_mass(gas);
ha1 = enthalpy_mass(gas);
Sa1 = sa1*mdot_air;
Xa1_int = mdot_air*exergy_mass(gas);         %%%%%%%%%%%%%%%%%%%%% a1 internal exergy
Xa1_flow = mdot_air*flowExergy_mass(gas);

% The second air point is at high pressure.
% Walk up in pressure adjusting via the polytropic efficiency.
%Pa2 = Pa1*Pressure_Ratio/Burner_Pressure_Ratio;
Pa2 = Pa1*Pressure_Ratio;
dP  = (Pa2 - Pa1)/steps;
s   = sa1;
h   = ha1;
for P = Pa1:dP:Pa2
    % Find the isentropic state at this pressure.
    set(gas,'S',s,'P',P);
    % Get the isentropic enthalpy difference (reversible work).
    hs = enthalpy_mass(gas);
    dhs = hs -h; 
    % The actual enthalpy difference is higher due to inefficiency.
    h = h + dhs/Air_Compressor_Polytropic_Efficiency;
    % Find the actual state at this pressure.
    set(gas,'H',h,'P',P);
    % Save the entropy for the next step in pressure.
    s = entropy_mass(gas);
end
Ta2 = temperature(gas);
ha2 = h;
sa2 = s;
Sa2 = sa2*mdot_air;
Xa2_int = mdot_air*exergy_mass(gas);  %%%%%%%%%%%%%%%%%%%%%%% a2 internal exergy
Xa2_flow = mdot_air*flowExergy_mass(gas); % a2 flow exergy


% Find the isentropic efficiency for this pressure ratio just for fun.
set(gas,'S',sa1,'P',Pa2);
ha2s = enthalpy_mass(gas);
Air_Compressor_Isentropic_Efficiency = (ha2s - ha1)/(ha2 - ha1);

% Start with fuel at ambient conditions.
xfuel = zeros(1,N);
xfuel(iCH4)  = 0.907;
xfuel(iC2H6) = 0.036;
xfuel(iC3H8) = 0.019;
xfuel(iN2)   = 0.018;
xfuel(iCO2)  = 0.010;
xfuel(iO2)   = 0.010;

% Find the lower heating value of the fuel.  Use extra oxygen to ensure
% complete combustion.  Since To is so low, dissociation is not an issue.
% Note that some versions of Cantera can't do 25C so use 300K.  Ug!
Nfuel = xfuel;                  % Use one mole of fuel
Noxid = 10;                     % Use excess oxygen (>2 for methane)
Nmix = Nfuel;                   % Put it in the mixture
Nmix(iO2) = Nmix(iO2) + Noxid;  % Add oxygen in excess of requirements
mass_mix = 0;                   % Find the mass of the mixture
for i=1:1:N
    mass_mix = mass_mix + Nmix(i)*M(i);
end
mass_fraction_O2 = Noxid*M(iO2)/mass_mix;   % Get the fraction of oxid
mass_fraction_fuel = 1 - mass_fraction_O2;  % Get the fraction of fuel

set(gas,'T',300,'P',101325,'X',Nmix);       % Put in the mixture
h_reactants = enthalpy_mass(gas);           % Find its enthalpy
equilibrate(gas,'TP');                      % Burn it at constant TP
h_products = enthalpy_mass(gas);            % Find the enthalpy
LHV_fuel = (h_reactants - h_products)/mass_fraction_fuel;    % Get LHV (mass)

Tf1 = To;
Pf1 = Po;
set(gas,'T',Tf1,'P',Pf1,'X',xfuel);
M_fuel = meanMolecularWeight(gas);
sf1 = entropy_mass(gas);
hf1 = enthalpy_mass(gas);
xf1_int = exergy_mass(gas);   %%%%%%%%%%%%%%%%% f1 specific internal exergy 
xf1_flow = flowExergy_mass(gas);

% The second fuel point is at high pressure (twice combustor pressure).
% Walk up in pressure adjusting via the polytropic efficiency.
Pf2 = 2*Pa2;
dP  = (Pf2 - Pf1)/steps;
s   = sf1;
h   = hf1;
for P = Pf1:dP:Pf2
    % Find the isentropic state at this pressure.
    set(gas,'S',s,'P',P);
    % Get the isentropic enthalpy difference (reversible work).
    hs = enthalpy_mass(gas);
    dhs = hs -h; 
    % The actual enthalpy difference is higher due to inefficiency.
    h = h + dhs/Fuel_Compressor_Polytropic_Efficiency;
    % Find the actual state at this pressure.
    set(gas,'H',h,'P',P);
    % Save the entropy for the next step in pressure.
    s = entropy_mass(gas);
end
Tf2 = temperature(gas);
hf2 = h;
sf2 = s;
xf2_int = exergy_mass(gas);  %%%%%%%%%%%%%%%% f2 specific internal exergy
xf2_flow = flowExergy_mass(gas);


% Find the isentropic efficiency for this pressure ratio just for fun.
set(gas,'S',sf1,'P',Pf2);
hf2s = enthalpy_mass(gas);
Fuel_Compressor_Isentropic_Efficiency = (hf2s - hf1)/(hf2 - hf1);

% Vary the fuel flow rate (at fixed air flow rate) to find the rate that
% gives the correct turbine inlet temperature.

% Choose a small step in fuel.
dfuel = (mdot_air/15)/2000;
for mdot_fuel=dfuel:dfuel:mdot_air
    mdot_mix = mdot_air + mdot_fuel;
    Sf1 = sf1*mdot_fuel;
    Sf2 = sf2*mdot_fuel;

    % Mix the air and fuel adiabatically.
    Pm3 = Pa2;
    hm3 = (mdot_fuel*hf2 + mdot_air*ha2)/mdot_mix;
    xmix = xfuel*mdot_fuel/M_fuel;
    xmix(iO2)  = xmix(iO2) + xair(iO2)*mdot_air/M_air;
    xmix(iN2)  = xmix(iN2) + xair(iN2)*mdot_air/M_air;           %%%% RH
    xmix(iAR)  = xair(iAR)*mdot_air/M_air;
    xmix(iCO2) = xair(iCO2)*mdot_air/M_air;
    xmix(iH2O) = xair(iH2O)*mdot_air/M_air;
%     xair(iN2) = (1-xair(iH2O))*0.757223/(1-0.031208);
%     xair(iO2) = (1-xair(iH2O))*0.202157/(1-0.031208);
%     xair(iAR) = (1-xair(iH2O))*0.009015/(1-0.031208);
%     xair(iCO2)= (1-xair(iH2O))*0.000397/(1-0.031208);
%     xair(iH2O)= 

    xmix = xmix/sum(xmix);
    set(gas,'H',hm3,'P',Pm3,'X',xmix);
    Tm3 = temperature(gas);
    sm3 = entropy_mass(gas);
    Sm3 = sm3*mdot_mix;
    Xm3_int = mdot_mix*exergy_mass(gas);       %%%%%%%%%%%%%%%% m3 internal exergy
    Xm3_flow = mdot_mix*flowExergy_mass(gas);  

    % Burn the mixture adiabatically.
    Pm4 = Pm3*Burner_Pressure_Ratio;
    set(gas,'H',hm3,'P',Pm4);
    equilibrate(gas,'HP');
    Tm4 = temperature(gas);
    
    % Look for a simple crossover since the steps are small and fast.
    if(Tm4 >= Tmax)
        Tpeak = Tm4;
        break
    end
end
sm4 = entropy_mass(gas);
hm4 = enthalpy_mass(gas);
Sm4 = sm4*mdot_mix;

xmix4(:) = moleFractions(gas);                                     %%%%%%%%%%%%%%%%%% moleFractions after combustion
Xm4_int = mdot_mix*exergy_mass(gas);            %%%%%%%%%%%%%%%%%%% m4 internal exergy
Xm4_flow = mdot_mix*flowExergy_mass(gas);

Xf1_int = mdot_fuel*xf1_int;                                          %%%%%%%%%%%%%%%%%% f1 f2 internal exergy
Xf2_int = mdot_fuel*xf2_int;
Xf1_flow = mdot_fuel*xf1_flow;                                          %%%%%%%%%%%%%%%%%% f1 f2 flow exergy
Xf2_flow = mdot_fuel*xf2_flow;


% Expand the products.  Use small steps and polytropic efficiency.
Pm5 = Po;
dP  = (Pm4 - Pm5)/steps;
s   = sm4;
h   = hm4;
for P = Pm4:-dP:Pm5
    % Find the starting state for the step.
    set(gas,'S',s,'P',P);
    % Find the isentropic enthapy change.
    hs = enthalpy_mass(gas);
    dhs = h - hs;
    % The actual change is less due to inefficiency.
    h = h - dhs*Mix_Turbine_Polytropic_Efficiency;
    % Find the actual state.
    set(gas,'H',h,'P',P);
    % Include the next line for shifting equilibrium.  Remove if you prefer
    % to leave the gas frozen in composition.  (A small difference.)
    equilibrate(gas,'HP');
    % Get the entropy for the next step.
    s = entropy_mass(gas);
end
Tm5 = temperature(gas);
hm5 = h;
sm5 = s;
Sm5 = sm5*mdot_mix;
xmix5=xmix4;
Xm5_int = mdot_mix*exergy_mass(gas);                %%%%%%%%%%%%%%%%%%%%% m5 internal exergy
Xm5_flow = mdot_mix*flowExergy_mass(gas);

% Find the isentropic efficiency for this pressure ratio just for fun.
set(gas,'S',sm4,'P',Pm5);
hm5s = enthalpy_mass(gas);
Mix_Turbine_Isentropic_Efficiency = (hm4 - hm5)/(hm4 - hm5s);

% Assemble an array of state points.
Tair  = [Ta1 Ta2];
Pair  = [Pa1 Pa2];
Sair  = [Sa1 Sa2];
Tfuel = [Tf1 Tf2];
Pfuel = [Pf1 Pf2];
Sfuel = [Sf1 Sf2];
Tmix  = [Tm4 Tm5];
Pmix  = [Pm4 Pm5];
Smix  = [Sm4 Sm5];

% Exergy List
x_int  = [Xf1_int Xf2_int Xa1_int Xa2_int Xm3_int Xm4_int Xm5_int]/mdot_fuel;      % semi-extensive internal exergy
x_flow = [Xf1_flow Xf2_flow Xa1_flow Xa2_flow Xm3_flow Xm4_flow Xm5_flow]/mdot_fuel;
X_int  = [Xf1_int Xf2_int Xa1_int Xa2_int Xm3_int Xm4_int Xm5_int];                % extensive exergy lists
X_flow = [Xf1_flow Xf2_flow Xa1_flow Xa2_flow Xm3_flow Xm4_flow Xm5_flow];

% Find the various energy transfers.
Air_Fuel_Mass_Ratio = mdot_air/mdot_fuel;
W_fuel_compressor   = mdot_fuel*(hf2 - hf1);
W_air_compressor    = mdot_air*(ha2 - ha1);
W_gross_gas_turbine = mdot_mix*(hm4 - hm5);
W_net_gas_turbine   = W_gross_gas_turbine - W_fuel_compressor - W_air_compressor;
Eff_gas_turbine     = W_net_gas_turbine/(mdot_fuel*LHV_fuel);                        % 1st law efficieny to LHV
Air_Specific_Work   = W_net_gas_turbine/mdot_air;

Eff_exergy        = (W_net_gas_turbine)/(Xf1_flow+Xa1_flow);             % Exergy efficiency

% Extensive exergy analysis
% X_des_total     = Xf1_flow+Xa1_flow-W_net_gas_turbine;
% X_des_fuel_comp = (Xf1_flow-Xf2_flow+W_fuel_compressor)/X_des_total*100;                   % Destroyed exergy at fuel compressor
% X_des_air_comp  = (Xa1_flow-Xa2_flow+W_air_compressor)/X_des_total*100;                    % Destroy exergy at air compressor
% X_des_premixer  = (Xf2_flow+Xa2_flow-Xm3_flow)/X_des_total*100;                            % Destroyed exergy at premixer
% X_des_combustor = (Xm3_flow-Xm4_flow)/X_des_total*100;                                     % Destroyed exergy at combustor
% X_des_turbine   = (Xm4_flow-Xm5_flow-W_gross_gas_turbine)/X_des_total*100;                 % Destroyed exergy at turbine
% X_des_exhaust   = Xm5_flow/X_des_total*100;
% X_des_list      = [X_des_fuel_comp X_des_air_comp X_des_premixer X_des_combustor X_des_turbine X_des_exhaust];

X_des_total     = Xf1_flow+Xa1_flow-W_net_gas_turbine;
X_des_fuel_comp = (Xf1_flow-Xf2_flow+W_fuel_compressor);                   % Destroyed exergy at fuel compressor
X_des_air_comp  = (Xa1_flow-Xa2_flow+W_air_compressor);                    % Destroy exergy at air compressor
X_des_premixer  = (Xf2_flow+Xa2_flow-Xm3_flow);                            % Destroyed exergy at premixer
X_des_combustor = (Xm3_flow-Xm4_flow);                                     % Destroyed exergy at combustor
X_des_turbine   = (Xm4_flow-Xm5_flow-W_gross_gas_turbine);                 % Destroyed exergy at turbine
X_des_exhaust   = Xm5_flow;
X_des_list      = [X_des_fuel_comp X_des_air_comp X_des_premixer X_des_combustor X_des_turbine X_des_exhaust];
X_des_list      = X_des_list/sum(X_des_list)*100;


% Exergy Destruction Plot
figure(1)
X_bar = categorical({'Fuel \n Compressor', 'Air \n Compressor', 'Premixer', 'Combustor', 'Turbine', 'Exhaust'});
X_bar = reordercats(X_bar,{'Fuel \n Compressor', 'Air \n Compressor', 'Premixer', 'Combustor', 'Turbine', 'Exhaust'});
bar(X_bar, X_des_list)
xlabel("Conversion Device")
ylabel("Exergy Loss [%]")
ylim([0 50])
text(1,30,['LHV Efficiency: ' string(Eff_gas_turbine*100)])
text(1,20,['Exergy Efficiency: ' string(Eff_exergy*100)])

figure(2)
station_list = [1 2 3 4 5]; % 1: Fuel-1 and air-1, 2: Fuel-2 and air-2, 3: Mixing (m3) 4: after combustion (m4), 5: after turbine (m5)
x_int_fuel = [x_int(1) x_int(2) x_int(5)]/10^6;
x_int_air = [x_int(3) x_int(4) x_int(5)]/10^6;
x_int_mix = [x_int(5) x_int(6) x_int(7)]/10^6;
plot(station_list(1:3), x_int_air, 'bo-')
hold on
plot(station_list(1:3), x_int_fuel, 'ro-')
hold on
plot(station_list(3:5), x_int_mix, 'ko-')
xlim([0.5 5.5])
ylim([0 70])

figure(3)
station_list = [1 2 3 4 5]; % 1: Fuel-1 and air-1, 2: Fuel-2 and air-2, 3: Mixing (m3) 4: after combustion (m4), 5: after turbine (m5)
x_flow_fuel = [x_flow(1) x_flow(2) x_flow(5)]/10^6;
x_flow_air = [x_flow(3) x_flow(4) x_flow(5)]/10^6;
x_flow_mix = [x_flow(5) x_flow(6) x_flow(7)]/10^6;

plot(station_list(1:3), x_flow_air, 'bo-')
hold on
plot(station_list(1:3), x_flow_fuel,'ro-')
hold on
plot(station_list(3:5), x_flow_mix, 'ko-')
xlim([0.5 5.5])
ylim([0 70])

%Make T-s plot.
figure(4)
plot(Sair/mdot_fuel/1000,Tair,'bo--')
hold on
plot([Sa2/mdot_fuel/1000 Sm3/mdot_fuel/1000],[Ta2 Tm3],'b--')
plot(Sfuel/mdot_fuel/1000,Tfuel,'ro--')
plot([Sf2/mdot_fuel/1000 Sm3/mdot_fuel/1000],[Tf2 Tm3],'r--')
plot(Smix/mdot_fuel/1000,Tmix,'ko--')
plot([Sm3/mdot_fuel/1000 Sm4/mdot_fuel/1000],[Tm3 Tm4],'ko--')
plot(Sm5/mdot_fuel/1000,Tm5,'ko')
text(Sa1/mdot_fuel/1000,Ta1,'a,1')
text(Sf1/mdot_fuel/1000,Tf1,'f,1')
text(Sa2/mdot_fuel/1000,Ta2,'a,2')
text(Sf2/mdot_fuel/1000,Tf2,'f,2')
text(Sm3/mdot_fuel/1000,Tm3,'m,3')
text(Sm4/mdot_fuel/1000,Tm4,'m,4')
text(Sm5/mdot_fuel/1000,Tm5,'m,5')
text(30,1750,sprintf('Thermal Efficiency:  %.1f%% (LHV)',100*Eff_gas_turbine))
text(30,1550,sprintf('Air-Specific Work:   %.2f MJ/kg-air',Air_Specific_Work/1e6))
text(30,1350,sprintf('Fuel LHV:                 %.1f MJ/kg',LHV_fuel/1e6))
scale = axis
axis([0 400 0 2000]);
xlabel('Fuel-Specific Entropy (kJ/kg_f_u_e_l-K)')
ylabel('Temperature (K)')
hold off
% plotfixer
legend('off')
