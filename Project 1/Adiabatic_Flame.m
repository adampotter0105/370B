% Calculate the adiabatic flame temperature and composition for isobaric or
% isochoric combustion of methane/air.
% C.F. Edwards, 12/28/13

clear
format compact
fprintf('*********************************************************\n')

% Make a Cantera gas object to work with.
gas = IdealGasMix('GRI30.xml')
nsp = nSpecies(gas);

% Set up some indices to the species of interest.
iCH4   = speciesIndex(gas,'CH4');
iO2    = speciesIndex(gas,'O2');
iN2    = speciesIndex(gas,'N2');

% Walk through a range of equivalence ratios.
% CH4:   fCH4 + 2(O2 + 3.76N2) -> fCO2 + 2fH2O + 2*3.76N2
for(i=1:1:101)
   f = 2*((i-1)/100)                % Set phi.
   phi(i) = f;
   x = zeros(nsp,1);                % Molefraction array with zeros.
   % Set up for methane.
   x(iCH4) = phi(i);                % Set fuel to phi.
   x(iO2)  = 2.0;                   % Set stoich oxygen.
   x(iN2)  = 2*3.76;                % And enough nitrogen to make air.
   % Set the composition, temperature, and pressure of reactants.
   set(gas,'Temperature',800,'Pressure',20*100000,'MoleFractions',x);
   % Choose one of the following lines to test your Cantera
   % installation.
%    equilibrate(gas,'HP');           % Burn them adiabatically.
   equilibrate(gas,'UV');           % Burn them adiabatically.
   tad(i) = temperature(gas);       % Save the temperature,
   xeq(:,i) = moleFractions(gas);   %   and the products too.
end

figure(1)
clf
plot(phi,tad);
xlabel('Equivalence Ratio');
ylabel('Adiabatic Flame Temperature (K)');
text(.1,2800,'800 K, 20 bar','FontSize',14);
title('Methane/Air')
plotfixer
legend('off')

figure(2)
clf
semilogy(phi,xeq);
axis([0 2 1.0e-3 1]);
% Use this loop to write a list that can help you sort the species
% Comment it out when you know what labels you want and where.
% for(i=1:1:nsp)
%     text(1,xeq(i,51),speciesName(gas,i),'FontSize',14);
% end
text(phi(50),xeq(speciesIndex(gas,'N2'),50),'N_2','FontSize',14);
text(phi(50),xeq(speciesIndex(gas,'CO2'),50),'CO_2','FontSize',14);
text(phi(60),xeq(speciesIndex(gas,'H2O'),60),'H_2O','FontSize',14);
text(phi(50),xeq(speciesIndex(gas,'CO'),50),'CO','FontSize',14);
text(phi(45),xeq(speciesIndex(gas,'OH'),45),'OH','FontSize',14);
text(phi(40),xeq(speciesIndex(gas,'NO'),40),'NO','FontSize',14);
text(phi(65),xeq(speciesIndex(gas,'H'),65),'H','FontSize',14);
text(phi(45),xeq(speciesIndex(gas,'O'),45),'O','FontSize',14);
text(phi(15),xeq(speciesIndex(gas,'O2'),15),'O_2','FontSize',14);
text(phi(60),xeq(speciesIndex(gas,'H2'),60),'H_2','FontSize',14);
xlabel('Equivalence Ratio');
ylabel('Mole Fraction');
title('Methane/Air')
plotfixer
legend('off')
