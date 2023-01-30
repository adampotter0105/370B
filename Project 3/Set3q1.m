% Set 3 Problem 1

clc;

T1 = 1600;
T2 = 1800;
T3 = 2000;
P_ratio = linspace(10,200);
gas = GRI30();

Eff_T1 = zeros(length(P_ratio),1);
Work_T1 = zeros(length(P_ratio),1);
for i = 1:length(P_ratio)
    [ef, wo] = gas_turbine_eff(gas, P_ratio(i), T1);
    Eff_T1(i,1) = ef;
    Work_T1(i,1) = wo;
end

Eff_T2 = zeros(length(P_ratio),1);
Work_T2 = zeros(length(P_ratio),1);
for i = 1:length(P_ratio)
    [ef, wo] = gas_turbine_eff(gas, P_ratio(i), T2);
    Eff_T2(i,1) = ef;
    Work_T2(i,1) = wo;
end

Eff_T3 = zeros(length(P_ratio),1);
Work_T3 = zeros(length(P_ratio),1);
for i = 1:length(P_ratio)
    [ef, wo] = gas_turbine_eff(gas, P_ratio(i), T3);
    Eff_T3(i,1) = ef;
    Work_T3(i,1) = wo;
end

P_ratio_scatter = linspace(10,200,5);
Eff_scat = zeros(length(P_ratio_scatter),3);
Work_scat = zeros(length(P_ratio_scatter),3);

for i = 1:3
    if i==1
        Temp = 1600;
    elseif i==2
        Temp = 1800;
    else
        Temp = 2000;
    end
    for j = 1:length(P_ratio_scatter)
        [ef, wo] = gas_turbine_eff(gas, P_ratio_scatter(j), Temp);
        Eff_scat(j,i) = ef;
        Work_scat(j,i) = wo;
    end
end


figure(1);
plot(Work_T1/1000,Eff_T1*100,'LineWidth',1.5,'Color','b');
hold on;
plot(Work_T2/1000,Eff_T2*100,'LineWidth',1.5,'Color','g');
plot(Work_T3/1000,Eff_T3*100,'LineWidth',1.5,'Color','r');

plot(Work_scat(:,1)/1000,Eff_scat(:,1)*100,'ob');
for i = 1:length(P_ratio_scatter)
    text(Work_scat(i,1)/1000,Eff_scat(i,1)*100,int2str(P_ratio_scatter(i)));
end
plot(Work_scat(:,2)/1000,Eff_scat(:,2)*100,'og');
for i = 1:length(P_ratio_scatter)
    text(Work_scat(i,2)/1000,Eff_scat(i,2)*100,int2str(P_ratio_scatter(i)));
end
plot(Work_scat(:,3)/1000,Eff_scat(:,3)*100,'or');
for i = 1:length(P_ratio_scatter)
    text(Work_scat(i,3)/1000,Eff_scat(i,3)*100,int2str(P_ratio_scatter(i)));
end

ylabel('First Law Efficiency (% LHV)');
xlabel('Air Specific Work (kJ/kg-air)');
title('Gas Turbine Carpet Plot');
legend('T = 1600 K', 'T = 1800 K', 'T = 2000 K');
hold off;



function [efficiency_lhv, air_specific_work] = gas_turbine_eff(gas, pressure_ratio, T_inlet)

    % Set the number of "differential" steps for the polytropic process.
    steps = 200;

    % Cycle specifications:
    To = 20+273.15;
    Po = 100000;
    Pressure_Ratio = pressure_ratio;
    Burner_Pressure_Ratio                 = 0.95;
    Air_Compressor_Polytropic_Efficiency  = 0.90;
    Fuel_Compressor_Polytropic_Efficiency = 0.90;
    Mix_Turbine_Polytropic_Efficiency     = 0.88;
    Tmax = T_inlet;

    % Make a gas object to work with in Cantera.
    % gas = Solution('gri30.xml');
    N = nSpecies(gas);
    iCH4 = speciesIndex(gas,'CH4');
    iC2H6 = speciesIndex(gas,'C2H6');
    iC3H8 = speciesIndex(gas,'C3H8');
    iCO2 = speciesIndex(gas,'CO2');
    iO2 = speciesIndex(gas,'O2');
    iN2 = speciesIndex(gas,'N2');
    M = molecularWeights(gas);

% Start with dry engineering air at ambient conditions.
    mdot_air = 1;
    xair = zeros(1,N);
    xair(iO2) = 0.21;
    xair(iN2) = 0.79;
    Ta1 = To;
    Pa1 = Po;
    set(gas,'T',Ta1,'P',Pa1,'X',xair);
    M_air = meanMolecularWeight(gas);
    sa1 = entropy_mass(gas);
    ha1 = enthalpy_mass(gas);
    Sa1 = sa1*mdot_air;

    % The second air point is at high pressure.
    % Walk up in pressure adjusting via the polytropic efficiency.
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

    % Find the isentropic efficiency for this pressure ratio just for fun.
%     set(gas,'S',sa1,'P',Pa2);
%     ha2s = enthalpy_mass(gas);
%     Air_Compressor_Isentropic_Efficiency = (ha2s - ha1)/(ha2 - ha1);

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

    % Find the isentropic efficiency for this pressure ratio just for fun.
%     set(gas,'S',sf1,'P',Pf2);
%     hf2s = enthalpy_mass(gas);
%     Fuel_Compressor_Isentropic_Efficiency = (hf2s - hf1)/(hf2 - hf1)

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
        xmix(iO2)  = xmix(iO2) + 0.21*mdot_air/M_air;
        xmix(iN2)  = xmix(iN2) + 0.79*mdot_air/M_air;
        xmix = xmix/sum(xmix);
        set(gas,'H',hm3,'P',Pm3,'X',xmix);
        Tm3 = temperature(gas);
        sm3 = entropy_mass(gas);
        Sm3 = sm3*mdot_mix;

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

    % Find the isentropic efficiency for this pressure ratio just for fun.
%     set(gas,'S',sm4,'P',Pm5);
%     hm5s = enthalpy_mass(gas);
%     Mix_Turbine_Isentropic_Efficiency = (hm4 - hm5)/(hm4 - hm5s);

    % Assemble an array of state points.
%     Tair  = [Ta1 Ta2];
%     Pair  = [Pa1 Pa2];
%     Sair  = [Sa1 Sa2];
%     Tfuel = [Tf1 Tf2];
%     Pfuel = [Pf1 Pf2];
%     Sfuel = [Sf1 Sf2];
%     Tmix  = [Tm4 Tm5];
%     Pmix  = [Pm4 Pm5];
%     Smix  = [Sm4 Sm5];

    % Find the various energy transfers.
    % Air_Fuel_Mass_Ratio = mdot_air/mdot_fuel;
    W_fuel_compressor   = mdot_fuel*(hf2 - hf1);
    W_air_compressor    = mdot_air*(ha2 - ha1);
    W_gross_gas_turbine = mdot_mix*(hm4 - hm5);
    W_net_gas_turbine   = W_gross_gas_turbine - W_fuel_compressor - W_air_compressor;
    efficiency_lhv     = W_net_gas_turbine/(mdot_fuel*LHV_fuel);
    air_specific_work   = W_net_gas_turbine/mdot_air;
end

