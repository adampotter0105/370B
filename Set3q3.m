% Set 3 Problem 3

clc;

gas = GRI30();

% Set the number of "differential" steps for the polytropic process.
steps = 200;

% Cycle specifications:
To = 25+273.15;
Po = 101325;
Pressure_Ratio = 12.2;
Burner_Pressure_Ratio                 = 0.95;
Air_Compressor_Polytropic_Efficiency  = 0.86;
Fuel_Compressor_Polytropic_Efficiency = 0.86;
Mix_Turbine_Polytropic_Efficiency     = 0.82;
Tmax = 1410;

HRSG_effectiveness = 0.90;
Econ_Pressure_Ratio = 0.92;
Boiler_Pressure_Ratio = 0.92;
Superheater_Pressure_Ratio = 0.92;
Condenser_Pressure_Ratio = 0.92;
Steam_Turbine_Polytropic_Efficiency = 0.75;
Condense_Pump_Polytropic_Efficiency = 0.85;
Feed_Pump_Polytropic_Efficiency = 0.85;

Pw6 = 6800;
Vaporfrac_6w = 0.88;
Pinch_Point_Temp_Diff = 20;
mdot_air = 144;

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
% mdot_air = 1;
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
Ha1 = ha1*mdot_air;

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
Ha2 = ha2*mdot_air;

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
    Hf1 = hf1*mdot_fuel;
    Hf2 = hf2*mdot_fuel;

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
    Hm3 = hm3*mdot_mix;

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
Hm4 = hm4*mdot_mix;

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
Hm5 = hm5*mdot_mix;
Pm5 = pressure(gas);

water = Water();

% initial cold storage state at ambient
Tw1 = To;
Pw1 = Po;
set(water,'T',Tw1,'P',Pw1);
hw1 = enthalpy_mass(water);
sw1 = entropy_mass(water);

Pw7 = Pw6*Condenser_Pressure_Ratio;

% known state at 6w
set(water, 'P', Pw6, 'Vapor', Vaporfrac_6w);
hw6 = enthalpy_mass(water);
Tw6 = temperature(water);
Tw7 = Tw6;
set(water, 'P', Pw7, 'T', Tw7);
hw7 = enthalpy_mass(water);
sw7 = entropy_mass(water);

Pw8 = Po;
[work78, hw8] = comp_work(Condense_Pump_Polytropic_Efficiency, Pw7, Pw8, hw7);
set(water, 'H', hw8, 'P', Pw8);
Tw8 = temperature(water);


Pw2_guess = 63000000;
Pw2 = Pw2_guess;

[work12, hw2] = comp_work(Feed_Pump_Polytropic_Efficiency, Pw1, Pw2, hw1);

set(water,'H',hw2,'P',Pw2);
Tw2 = temperature(water);
sw2 = enthalpy_mass(water);

Pw3 = Econ_Pressure_Ratio*Pw2;
Pw4 = Boiler_Pressure_Ratio*Pw3;
Pw5 = Superheater_Pressure_Ratio*Pw4;

Pm6 = Pm5;
Pm7 = Pm6;
Pm8 = Pm7;

% Find heat transferred through HRSG
set(gas, 'T',Tm5,'P',Pw5);
hlc1 = enthalpy_mass(gas);
set(gas,'T',Tw2,'P',Pw2);
hlc2 = enthalpy_mass(gas);
qmax = hlc1 - hlc2;
qactual = HRSG_effectiveness*qmax;

hm8 = hm5 - qactual;
set(gas,'H',hm8,'P',Pm8);
Tm8 = temperature(gas);

hw5 = hw2 + qactual;
set(water, 'H',hw5,'P',Pw5);
Tw5 = temperature(water);

% Pw5_guess = 49000000;
% Tw5 = 750;
% set(water,'T',Tw5,'P',Pw5_guess);
% hw5 = enthalpy_mass(water);

[work56, hw6hold, X6] = turb_work(Steam_Turbine_Polytropic_Efficiency, Pw5, Pw6, hw5);
set(water, 'H',hw6hold,'P',Pw6);
disp(X6);

offset = abs(Hm5);
Ha1_off = Ha1 + offset;
Ha2_off = Ha2 + offset;
Hf1_off = Hf1 + offset;
Hf2_off = Hf2 + offset;
Hm3_off = Hm3 + offset;
Hm4_off = Hm4 + offset;
Hm5_off = Hm5 + offset;

% Find the isentropic efficiency for this pressure ratio just for fun.
%     set(gas,'S',sm4,'P',Pm5);
%     hm5s = enthalpy_mass(gas);
%     Mix_Turbine_Isentropic_Efficiency = (hm4 - hm5)/(hm4 - hm5s);

% Assemble an array of state points.
Tair  = [Ta1 Ta2];
Pair  = [Pa1 Pa2];
Sair  = [Sa1 Sa2];
Hair_off  = [Ha1_off Ha2_off];
Tfuel = [Tf1 Tf2];
Pfuel = [Pf1 Pf2];
Sfuel = [Sf1 Sf2];
Hfuel_off = [Hf1_off Hf2_off];
Tmix  = [Tm4 Tm5];
Pmix  = [Pm4 Pm5];
Smix  = [Sm4 Sm5];
Hmix_off  = [Hm4_off Hm5_off];

% Find the various energy transfers.
% Air_Fuel_Mass_Ratio = mdot_air/mdot_fuel;
W_fuel_compressor   = mdot_fuel*(hf2 - hf1);
W_air_compressor    = mdot_air*(ha2 - ha1);
W_gross_gas_turbine = mdot_mix*(hm4 - hm5);
W_net_gas_turbine   = W_gross_gas_turbine - W_fuel_compressor - W_air_compressor;
efficiency_lhv     = W_net_gas_turbine/(mdot_fuel*LHV_fuel);
air_specific_work   = W_net_gas_turbine/mdot_air;

figure(1)
plot(Hair_off/mdot_fuel/1e6,Tair,'bo--')
hold on
plot([Ha2_off/mdot_fuel/1e6 Hm3_off/mdot_fuel/1e6],[Ta2 Tm3],'b--')
plot(Hfuel_off/mdot_fuel/1e6,Tfuel,'ro--')
plot([Hf2_off/mdot_fuel/1e6 Hm3_off/mdot_fuel/1e6],[Tf2 Tm3],'r--')
plot(Hmix_off/mdot_fuel/1e6,Tmix,'ko--')
plot([Hm3_off/mdot_fuel/1e6 Hm4_off/mdot_fuel/1e6],[Tm3 Tm4],'ko--')
plot(Hm5_off/mdot_fuel/1e6,Tm5,'ko')
text(Ha1_off/mdot_fuel/1e6,Ta1,'a,1')
text(Hf1_off/mdot_fuel/1e6,Tf1,'f,1')
text(Ha2_off/mdot_fuel/1e6,Ta2,'a,2')
text(Hf2_off/mdot_fuel/1e6,Tf2,'f,2')
text(Hm3_off/mdot_fuel/1e6,Tm3,'m,3')
text(Hm4_off/mdot_fuel/1e6,Tm4,'m,4')
text(Hm5_off/mdot_fuel/1e6,Tm5,'m,5')
% text(30,1750,sprintf('Thermal Efficiency:  %.1f%% (LHV)',100*Eff_gas_turbine))
% text(30,1550,sprintf('Air-Specific Work:   %.2f MJ/kg-air',Air_Specific_Work/1e6))
% text(30,1350,sprintf('Fuel LHV:                 %.1f MJ/kg',LHV_fuel/1e6))

xlabel('Fuel-Specific Enthalpy Difference (MJ/kg_f_u_e_l-K)')
ylabel('Temperature (K)')

function [total_w, H] = comp_work(eta, Pin, Pout, h0)
    w = Water();
    set(w, "H", h0, "P", Pin) % Initial Conditions of fluid
    n = 101; % discretization points for pressure
    h = (Pout-Pin)/n; % integration step
    total_w = 0; % total work output of turbine
    for p = Pin+h:h:Pout
        h1 = enthalpy_mass(w);
        s0 = entropy_mass(w);
        
        setState_SP(w, [s0, p]) % isentropic expansion
        h2 = enthalpy_mass(w);
        dw_s = h2-h1; % isentropic work
    
        dw = dw_s/eta; % actual work
        setState_HP(w, [h1+dw, p]) % actual new state of steam 
        total_w = total_w + dw; % sum all steps
    end
    H = enthalpy_mass(w);
end

function [total_w, H, X] = turb_work(eta, Pin, Pout, h0)
    w = Water;
    set(w, "H", h0, "P", Pin) % Initial Conditions of water
    n = 100; % discretization points for pressure
    h = (Pout-Pin)/n; % integration step
    total_w = 0; % total work output of turbine
    for p = Pin+h:h:Pout
        h1 = enthalpy_mass(w);
        s0 = entropy_mass(w);
        setState_SP(w, [s0, p]) % isentropic expansion
        h2 = enthalpy_mass(w);
        dw_s = h2-h1; % isentropic work
    
        dw = dw_s*eta; % actual work
        setState_HP(w, [h1+dw, p]) % actual new state of steam 
        total_w = total_w + dw; % sum all steps
    end
    X = vaporFraction(w);
    H = enthalpy_mass(w);
end


