% ME370B Project 4 - Part 1
% Winter 2023
% Andy Huynh

clear

global T0 P0 nOCR
T0 = 450+273.15; % pre-heated T
P0 = 10*100000; % pre-pressurized P

nOCR = 100; % number of O/C ratio points
OCR = linspace(0,1.0,nOCR); % O/C ratio range

%% Calculate Directly for Methane
gas = Solution('gasification_small.xml'); % set mechanism
nsp = nSpecies(gas); % number of species in mechanism
species_names = speciesNames(gas); % get list of species names

% Find species indices
ich4 = speciesIndex(gas,'CH4');
ih2o = speciesIndex(gas,'H2O');
io2 = speciesIndex(gas,'O2');
in2 = speciesIndex(gas,'N2');
ico = speciesIndex(gas,'CO');
ih2 = speciesIndex(gas,'H2');

Xin_direct(nsp,nOCR) = 0; % initialize matrix to store initial mole fractions
Teq_direct(nOCR) = 0; % initialize vector to store temperature at equilibrium
Xeq_direct(nsp,nOCR) = 0; % intialize matrix to store mole fractions at equilibrium

for i = 1:nOCR
    x = zeros(nsp,1); % intialize matrix to store mole fraction data
    x(ich4,1) = 1.0; % set mole fractions
    x(ih2o,1) = 1.0;
    x(io2,1) = OCR(i);
    x(in2,1) = 3.76*OCR(i);
    Xin_direct(:,i) = x;
    set(gas,'T',T0,'P',P0,'X',x); % define gas properties
    equilibrate(gas,'HP'); % autothermal reformer
    Teq_direct(i) = temperature(gas); % temperature at equilibrium (K)
    Xeq_direct(:,i) = moleFractions(gas); % mole fractions at equilibrium
end

%% CALCULATE INDIRECT
% Set target fuel values
% [methane, methanol, gasoline]
HC_fuels = [4 4 1.95]; % TODO: Find real values !!!!!!!!
OC_fuels = [0 1/6 0];
LHV_fuels = [50 19.9 43.4]; % MJ/kg  (engineering toolbox)
Cp_fuels = [3.6 3.62 3.6]; % 
nFuels = max(size(HC_fuels));

% Initilize Data Logging Variables
Xin_indirect(nFuels, nsp, nOCR) = 0;
Xeq_indirect(nFuels, nsp, nOCR) = 0;
Teq_indirect(nFuels, nOCR) = 0;

for f = 1:nFuels
    for ocr = 1:nOCR
        [X_in, X_out, T_out] = SMR(HC_fuels(f), OC_fuels(f), OCR(ocr), LHV_fuels(f), Cp_fuels(f));
        Xin_indirect(f,:,ocr) = X_in;
        Xeq_indirect(f,:,ocr) = X_out;
        Teq_indirect(f,ocr)= T_out;
    end
end

%% GENERATE PLOTS

% Compare direct and indirect: Plot temperature vs. O/C molar feed ratio
figure(1)
hold on
plot(OCR,Teq_direct-273.15);
plot(OCR, Teq_indirect(1,:)-273.15);
xlabel('Oxygen/Carbon Molar Feed Ratio');
ylabel('Equilibrium Temperature (°C)');
legend(["Direct CH4", "Indirect CH4"])
xticks([0 0.2 0.4 0.6 0.8 1.0]);
hold off

% Compare fuels: Plot temperature vs. O/C molar feed ratio
figure(2)
hold on
plot(OCR,Teq_indirect(1,:)-273.15);
plot(OCR,Teq_indirect(2,:)-273.15);
plot(OCR,Teq_indirect(3,:)-273.15);
xlabel('Oxygen/Carbon Molar Feed Ratio');
ylabel('Equilibrium Temperature (°C)');
legend(["Methane", "Methanol", "Gasoline"])
xticks([0 0.2 0.4 0.6 0.8 1.0]);
hold off

% Only plot species with over 0.1% mole fraction
[Xeq_direct_maj, species_direct] = majSpecies(Xeq_direct, species_names);
[Xeq_indirect_maj1, species_indirect1] = majSpecies(squeeze(Xeq_indirect(1,:,:)), species_names);
[Xeq_indirect_maj2, species_indirect2] = majSpecies(squeeze(Xeq_indirect(2,:,:)), species_names);
[Xeq_indirect_maj3, species_indirect3] = majSpecies(squeeze(Xeq_indirect(3,:,:)), species_names);

% TODO: FIX LEGENDS!!!
% Plot major species vs. O/C molar feed ratio
figure(3)
hold on
plot(OCR,Xeq_direct_maj, "-");
plot(OCR,Xeq_indirect_maj1, "--");
xlabel('Oxygen/Carbon Molar Feed Ratio');
ylabel('Equilibrium Mole Fraction');
legend(["Direct CH4", "Indirect CH4"])
xticks([0 0.2 0.4 0.6 0.8 1.0]);
legend(species_direct)
hold off

figure(4)
hold on
plot(OCR,Xeq_indirect_maj1, "-");
plot(OCR,Xeq_indirect_maj2, "--");
plot(OCR,Xeq_indirect_maj3, ".");
xlabel('Oxygen/Carbon Molar Feed Ratio');
ylabel('Equilibrium Mole Fraction');
legend(["Methane", "Methanol", "Gasoline"])
xticks([0 0.2 0.4 0.6 0.8 1.0]);
legend(species_indirect3)
hold off

function [Xeq_major, species_names_maj] = majSpecies(X_in, species_names)
    nsp = size(X_in,1);
    global nOCR
    % Track major species with more than 0.1% mole fraction
    Xeq_major = X_in; % store major species only
    species_names_maj = species_names;
    for n = nsp:-1:1
        if Xeq_major(n,nOCR/2) < 0.001 % get equilibrium composition at midpoint
            Xeq_major(n,:) = []; % delete rows
            species_names_maj(n) = [];
        end
    end
end

function [X_in, X_out, T_out] = SMR(HC_fuel, OC_fuel, OC_feed, LHV, Cp)
    global T0 P0
    gas = Solution('gasification_small.xml');
    nsp = nSpecies(gas); % number of species in gas
    x = zeros(nsp,1); % intialize matrix to store mole fraction data

    % Find species indices
    ico = speciesIndex(gas,'CO');
    ih2 = speciesIndex(gas,'H2');
    io2 = speciesIndex(gas,'O2');
    in2 = speciesIndex(gas,'N2');
    ih2o = speciesIndex(gas, 'H2O');
    ico2 = speciesIndex(gas, 'CO2');
    
    % C-H2O ratio assumed to be 1
    x(ico,1) = 1.0; % C with O from water
    x(ih2,1) = HC_fuel/2 + 1; % including H from water
    x(io2,1) = (OC_fuel + OC_feed)/2; % O from fuel and air
    x(in2,1) = (OC_feed)*3.76/2;
    X_in = x;
    mass_tot = [28 2 16 14]*[x(ico,1) x(ih2,1) x(io2,1) x(in2,1)]';
    mass_frac_fuel = (  (1)*12 + (HC_fuel/2)*2 +  (OC_fuel/2)*16  )/mass_tot;
    mass_frac_air_water = ( (OC_feed/2)*16 + (OC_feed*3.76/2)*14 + 1*18 )/mass_tot;

    % Setting up gas with raw elements
    mass_fuel = 12*1 + 1*HC_fuel + 8*OC_fuel;
    % Enthalpy contribution from fuel
    x_comb = zeros(nsp,1); 
    x_comb(ico2,1) = 1;
    x_comb(ih2o,1) = HC_fuel/2;
    mass_frac_comb = (44*x_comb(ico2,1) + 18*x_comb(ih2o,1)  )/mass_fuel;
    set(gas, 'T', 298, 'P', oneatm, 'X', x_comb);
    h_comb_products = enthalpy_mass(gas); 
    % Enthalpy of oxygen input
    x_o2 = zeros(nsp,1); 
    x_o2(io2,1) = 1;
    mass_frac_o2 = (2*1 + 0.5*HC_fuel - OC_fuel)*16/mass_fuel; % o2 needed to combust 1 carbon worth of fuel
    set(gas, 'T', 298, 'P', oneatm, 'X', x_o2); % TODO: finish converting LHV to enthalpy
    h_o2 = enthalpy_mass(gas);
    h_fuel = LHV*1e6 + mass_frac_comb*h_comb_products - mass_frac_o2*h_o2; % fuel + air = LHV + products

    % Get enthalpy contribution from water and air
    x_aw = zeros(nsp,1); 
    x_aw(io2,1) = OC_feed/2;
    x_aw(in2,1) = (OC_feed)*3.76/2;
    x_aw(ih2o) = 1;
    set(gas, 'T', T0, 'P', P0, 'X', x_aw);
    h_air_water = enthalpy_mass(gas); 
    
    h_in = mass_frac_fuel*h_fuel + mass_frac_air_water*h_air_water; % h_fuel at STP
    %actual enthalpy_mass of methane: -4.6493e+06

    set(gas,'H',h_in,'P',P0,'X',X_in)
    equilibrate(gas,'HP'); % autothermal reformer
    T_out = temperature(gas); % temperature at equilibrium (K)
    X_out = moleFractions(gas); % mole fractions at equilibrium

end