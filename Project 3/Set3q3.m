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
Xw7 = 0; % saturated liquid comes out of condenser
set(water, 'P', Pw7, 'Vapor', Xw7);
hw7 = enthalpy_mass(water);
sw7 = entropy_mass(water);
Tw7 = temperature(water);

Pw8 = Po;
[work78, hw8] = comp_work(Condense_Pump_Polytropic_Efficiency, Pw7, Pw8, hw7);
set(water, 'H', hw8, 'P', Pw8);
Tw8 = temperature(water);

% took out iteration loop to save run time here; this pressure was found
% through iteration
Pw2_guess = 10490000;
Pw2 = Pw2_guess;

[work12, hw2] = comp_work(Feed_Pump_Polytropic_Efficiency, Pw1, Pw2, hw1);

set(water,'H',hw2,'P',Pw2);
Tw2 = temperature(water);
sw2 = enthalpy_mass(water);

Pw3 = Econ_Pressure_Ratio*Pw2;
Pw4 = Boiler_Pressure_Ratio*Pw3;
Pw5 = Superheater_Pressure_Ratio*Pw4;

VapFracw3 = 0;
set(water,'P',Pw3,'Vapor',VapFracw3);
hw3 = enthalpy_mass(water);
Tw3 = temperature(water);

VapFracw4 = 1;
set(water,'T',Tw3,'Vapor',VapFracw4);
hw4 = enthalpy_mass(water);
Tw4 = temperature(water);

Pm6 = Pm5;
Pm7 = Pm6;
Pm8 = Pm7;

% Find heat transferred through HRSG
set(water, 'T',Tm5,'P',Pw5);
hlc1 = enthalpy_mass(water);
set(water,'T',Tw2,'P',Pw2);
hlc2 = enthalpy_mass(water);
qmax = hlc1 - hlc2;
qactual = HRSG_effectiveness*qmax;

% hm8 = hm5 - qactual;
% set(gas,'H',hm8,'P',Pm8,'X');
% Tm8 = temperature(gas);

hw5 = hw2 + qactual;
set(water, 'H',hw5,'P',Pw5);
Tw5 = temperature(water);

[work56, hw6hold, X6] = turb_work(Steam_Turbine_Polytropic_Efficiency, Pw5, Pw6, hw5);
set(water, 'H',hw6hold,'P',Pw6);
% disp(X6);

% find mass flow rate through pinch point analysis
mdot_water_guess = 20;
mdot_water = mdot_water_guess;
while true
    
    m57diff = (mdot_water*(hw5-hw3))/mdot_mix;
    set(gas, 'P',Pm7,'H',hm5-m57diff);
    Tm7 = temperature(gas);
    pinch = Tm7-Tw3;
    if pinch - Pinch_Point_Temp_Diff <= 0.0001
        break
    end
    mdot_water = mdot_water + 0.001;
end

m56diff = mdot_water*(hw5-hw4)/mdot_mix;
hm6 = hm5 - m56diff;
set(gas,'P',Pm6,'H',hm6);
Tm6 = temperature(gas);
Hm6 = hm6*mdot_mix;

hm7 = hm5 - m57diff;
set(gas,'P',Pm7,'H',hm7);
Tm7 = temperature(gas);
Hm7 = hm7*mdot_mix;

m78diff = mdot_water*(hw3-hw2)/mdot_mix;
hm8 = hm7 - m78diff;
set(gas,'P',Pm8,'H',hm8);
Tm8 = temperature(gas);
Hm8 = hm8*mdot_mix;

Hw1 = mdot_water*hw1;
Hw2 = mdot_water*hw2;
Hw3 = mdot_water*hw3;
Hw4 = mdot_water*hw4;
Hw5 = mdot_water*hw5;
Hw6 = mdot_water*hw6;
Hw7 = mdot_water*hw7;
Hw8 = mdot_water*hw8;

% Calculate Power produced by each turbine
W_gas_turb = mdot_mix*(hm4-hm5);
W_fuel_comp = mdot_fuel*(hf2 - hf1);
W_air_comp = mdot_air*(ha2 - ha1);
W_gas_net = W_gas_turb - W_fuel_comp - W_air_comp;

W_steam_turb = mdot_water*(hw5-hw6);
W_feed_pump = mdot_water*(hw2-hw1);
W_cond_pump = mdot_water*(hw8-hw7);
W_steam_net = W_steam_turb - W_feed_pump - W_cond_pump;
    
offset = abs(Hm5);
Ha1_off = Ha1 + offset;
Ha2_off = Ha2 + offset;
Hf1_off = Hf1 + offset;
Hf2_off = Hf2 + offset;
Hm3_off = Hm3 + offset;
Hm4_off = Hm4 + offset;
Hm5_off = Hm5 + offset;
Hm6_off = Hm6 + offset;
Hm7_off = Hm7 + offset;
Hm8_off = Hm8 + offset;

off_wat = abs(Hw5);
Hw1_off = Hw1 + off_wat;
Hw2_off = Hw2 + off_wat;
Hw3_off = Hw3 + off_wat;
Hw4_off = Hw4 + off_wat;
Hw5_off = Hw5 + off_wat;
Hw6_off = Hw6 + off_wat;
Hw7_off = Hw7 + off_wat;
Hw8_off = Hw8 + off_wat;

% Assemble an array of state points.
Tair  = [Ta1 Ta2];
Pair  = [Pa1 Pa2];
Sair  = [Sa1 Sa2];
Hair_off  = [Ha1_off Ha2_off];
Tfuel = [Tf1 Tf2];
Pfuel = [Pf1 Pf2];
Sfuel = [Sf1 Sf2];
Hfuel_off = [Hf1_off Hf2_off];
Tmix  = [Tm4 Tm5 Tm6 Tm7 Tm8];
% Pmix  = [Pm4 Pm5];
% Smix  = [Sm4 Sm5];
Hmix_off  = [Hm4_off Hm5_off Hm6_off Hm7_off Hm8_off];

Twat = [Tw1 Tw2 Tw3 Tw4 Tw5 Tw6 Tw7 Tw8];
Hwat_off = [Hw1_off Hw2_off Hw3_off Hw4_off Hw5_off Hw6_off Hw7_off Hw8_off];


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
text(Hm6_off/mdot_fuel/1e6,Tm6,'m,6')
text(Hm7_off/mdot_fuel/1e6,Tm7,'m,7')
text(Hm8_off/mdot_fuel/1e6,Tm8,'m,8')

plot(Hwat_off/mdot_fuel/1e6,Twat,'mo--');
text(Hw1_off/mdot_fuel/1e6,Tw1,'w,1');
text(Hw2_off/mdot_fuel/1e6,Tw2,'w,2');
text(Hw3_off/mdot_fuel/1e6,Tw3,'w,3');
text(Hw4_off/mdot_fuel/1e6,Tw4,'w,4');
text(Hw5_off/mdot_fuel/1e6,Tw5,'w,5');
text(Hw6_off/mdot_fuel/1e6,Tw6,'w,6');
text(Hw7_off/mdot_fuel/1e6,Tw7,'w,7');
text(Hw8_off/mdot_fuel/1e6,Tw8,'w,8');

Wg = round(W_gas_net/1e6,4,'significant');
Ws = round(W_steam_net/1e6,4,'significant');
text(-28,1500,strcat('Gas Turbine Power:  ', num2str(Wg), ' MW'),'FontWeight','bold');
text(-28,1400,strcat('Steam Turbine Power:  ', num2str(Ws), ' MW'),'FontWeight','bold');

xlabel('Fuel-Specific Enthalpy Difference (MJ/kg_f_u_e_l-K)');
ylabel('Temperature (K)');
title('Combined Cycle');

% add vapor dome
% Set up a fluid object to work with in Cantera.
fluid = Water();

% Get the critical point props.
Tc = critTemperature(water);
Pc = critPressure(water);
set(water,'T',Tc,'P',Pc);
sc = entropy_mass(water);
uc = intEnergy_mass(water);
hc = enthalpy_mass(water);
% Set the triple point props.  Use a value epsilon above Tt to get the
% bottom edge of the vapor dome.
Tt = 273.165;
Pt = satPressure(water,Tt);
setState_Tsat(water,[Tt 0]);
sft = entropy_mass(water);
uft = intEnergy_mass(water);
hft = enthalpy_mass(water);
vft = 1/density(water);
setState_Tsat(water,[Tt 1]);
sgt = entropy_mass(water);
ugt = intEnergy_mass(water);
hgt = enthalpy_mass(water);
vgt = 1/density(water);
% Set the limits for data curves.
Tmin = Tt;
Tmax = 1100;
Pmin = Pt;
Pmax = 10000*1e5;
% Make a vapor dome.
dT = (Tc-Tmin)/150;
i = 1;
%setState_Tsat(water,[Tmin 0]);
for T=Tmin:dT:Tc-dT    % Stop short of the critical point.
    Tsatline(i) = T;
    setState_Tsat(water,[T 0]);
    sliqline(i) = entropy_mass(water);
    uliqline(i) = intEnergy_mass(water);
    hliqline(i) = enthalpy_mass(water);
    setState_Tsat(water,[T 1]);
    svapline(i) = entropy_mass(water);
    uvapline(i) = intEnergy_mass(water);
    hvapline(i) = enthalpy_mass(water);
    i = i+1;
end
Tsatline(i) = Tc;   % Add the critical point now.
sliqline(i) = sc;
uliqline(i) = uc;
hliqline(i) = hc;
svapline(i) = sc;
uvapline(i) = uc;
hvapline(i) = hc;

hliqline = (hliqline*mdot_water + off_wat)/mdot_fuel;
hvapline = (hvapline*mdot_water + off_wat)/mdot_fuel;
plot(hliqline/1e6,Tsatline,'k')
plot(hvapline/1e6,Tsatline,'k')
hold off;

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


