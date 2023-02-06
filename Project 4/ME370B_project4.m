% ME370B Project 4 - Part 1
% Winter 2023
% Andy Huynh

gas = Solution('gasification_small.xml'); % set mechanism
nsp = nSpecies(gas); % number of species in mechanism

T0 = 450+273.15; % pre-heated T
P0 = 10*100000; % pre-pressurized P

% Find species indices
ich4 = speciesIndex(gas,'CH4');
ih2o = speciesIndex(gas,'H2O');
io2 = speciesIndex(gas,'O2');
in2 = speciesIndex(gas,'N2');
ico = speciesIndex(gas,'CO');
ih2 = speciesIndex(gas,'H2');

nOCR = 100; % number of O/C ratio points
OCR = linspace(0,1.0,nOCR); % O/C ratio range
Xin(nsp,nOCR) = 0; % initialize matrix to store initial mole fractions
Teq(nOCR) = 0; % initialize vector to store temperature at equilibrium
Xeq(nsp,nOCR) = 0; % intialize matrix to store mole fractions at equilibrium
LHVin = 0; % initialize vector to store LHV before reforming
LHVout = 0; % initilize vector to store LHV after reforming
Yin(nsp,nOCR) = 0; % initialize matrix to store mass fractions before reforming
Yout(nsp,nOCR) = 0; % initialize matrix to store mass fractions after reforming
EXin = 0; % initialize vector to store exergy going in
EXout = 0; % initialize vector to store exergy going in
MWin = 0; % intialize vector to store molecular weights going in
MWout = 0; % initialize vector to store molecular weights going out


for i = 1:nOCR
    x = zeros(nsp,1); % intialize matrix to store mole fraction data
    x(ich4,1) = 1.0; % set mole fractions
    x(ih2o,1) = 1.0;
    x(io2,1) = OCR(i);
    x(in2,1) = 3.76*OCR(i);
    Xin(:,i) = x;
    set(gas,'T',T0,'P',P0,'X',x); % define gas properties
    LHVin(i) = LHV_mass(gas); % LHV going in (J/kg)
    Yin(:,i) = massFractions(gas); % mass fractions going in
    %EXin(i) = flowExergy_mass(gas); % exergy going in
    MWin(i) = meanMolecularWeight(gas); % MW going in
    equilibrate(gas,'HP'); % autothermal reformer
    Teq(i) = temperature(gas); % temperature at equilibrium (K)
    Xeq(:,i) = moleFractions(gas); % mole fractions at equilibrium
    LHVout(i) = LHV_mass(gas); % LHV going out (J/kg)
    Yout(:,i) = massFractions(gas); % mass fractions going out
    %EXout(i) = flowExergy_mass(gas); % exergy going out
    MWout(i) = meanMolecularWeight(gas); % MW going out
end

% Plot temperature vs. O/C molar feed ratio
figure(1)
hold on
plot(OCR,Teq-273.15);
xlabel('Oxygen/Carbon Molar Feed Ratio');
ylabel('Equilibrium Temperature (°C)');
xticks([0 0.2 0.4 0.6 0.8 1.0]);
yticks([200 400 600 800 1000 1200 1400 1600]);
hold off

% Track major species with more than 0.1% mole fraction
for i = 1:nsp
    species_names(i) = speciesName(gas,i); % get list of species names
end

Xeq_major = Xeq; % store major species only
for n = nsp:-1:1
    if Xeq_major(n,nOCR/2) < 0.001 % get equilibrium composition at midpoint
        Xeq_major(n,:) = []; % delete rows
        species_names(n) = [];
    end
end

% Plot major species vs. O/C molar feed ratio
figure(2)
hold on
plot(OCR,Xeq_major);
xlabel('Oxygen/Carbon Molar Feed Ratio');
ylabel('Equilibrium Mole Fraction');
xticks([0 0.2 0.4 0.6 0.8 1.0]);
yticks([0 0.1 0.2 0.3 0.4 0.5]);
legend(species_names)
hold off

% Calculate cold-gas efficiency, exergy efficiency, and syngas yield
for i = 1:nOCR
    CGE(i) = LHVout(i)/LHVin(i)*100;
    %XE(i) = EXout(i)/EXin(i)*100;
    SGY(i) = (Xeq(ico,i)+Xeq(ih2,i))/Xin(ich4,i)*sum(Xin(:,i))*MWin(i)/MWout(i);
end

% Plot CGE,XE, and SGY vs. O/C molar feed ratio
figure(3)
hold on
plot(OCR,CGE,'DisplayName','Cold Gas (LHV)');
%plot(OCR,XR,'DisplayName','Exergy');
xlabel('Oxygen/Carbon Molar Feed Ratio');
ylabel('Efficiency (%)');
xticks([0 0.2 0.4 0.6 0.8 1.0]);
%yticks([60 70 80 90 100 110]);
legend()
hold off

figure(4)
hold on
plot(OCR,SGY);
xlabel('Oxygen/Carbon Molar Feed Ratio');
ylabel('Syngas Yield');
xticks([0 0.2 0.4 0.6 0.8 1.0]);
yticks([0 0.5 1.0 1.5 2.0 2.5 3.0]);
hold off