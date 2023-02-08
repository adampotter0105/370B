% Set 4 Problem 1

clc;

global mu_o To Po xo

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
ico2 = speciesIndex(gas,'CO2');
iar = speciesIndex(gas,'AR');

% dead state
To = 298.15;
Po = 101325;
xo = zeros(1,nsp);
xo(in2)  = 0.757223;
xo(io2)  = 0.202157;
xo(ih2o) = 0.031208;
xo(iar)  = 0.009015;
xo(ico2) = 0.000397;
set(gas,'T',To,'P',Po,'X',xo);
mu_o = chemPotentials(gas);

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
    resc = norm(x);
    x = x/resc;
    Xin(:,i) = x;
    set(gas,'T',T0,'P',P0,'X',x); % define gas properties
    LHVin(i) = LHV_mass(gas); % LHV going in (J/kg)
    Yin(:,i) = massFractions(gas); % mass fractions going in
    EXin(i) = flowExergy_mass(gas); % exergy going in
    
    MWin(i) = meanMolecularWeight(gas); % MW going in
    equilibrate(gas,'HP'); % autothermal reformer
    Teq(i) = temperature(gas); % temperature at equilibrium (K)
    Xeq(:,i) = moleFractions(gas); % mole fractions at equilibrium
    LHVout(i) = LHV_mass(gas); % LHV going out (J/kg)
    Yout(:,i) = massFractions(gas); % mass fractions going out
    % EXout(i) = flowExergy_mass(gas); % exergy going out
    MWout(i) = meanMolecularWeight(gas); % MW going out
    EXout(i) = flowExergy_mass(gas);
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
    XE(i) = EXout(i)/EXin(i)*100;
    SGY(i) = (Xeq(ico,i)+Xeq(ih2,i))/Xin(ich4,i)*sum(Xin(:,i))*MWin(i)/MWout(i);
end

% Plot CGE,XE, and SGY vs. O/C molar feed ratio
figure(3)
hold on
plot(OCR,CGE,'DisplayName','Cold Gas (LHV)');
plot(OCR,XE,'DisplayName','Exergy');
xlabel('Oxygen/Carbon Molar Feed Ratio');
ylabel('Efficiency (%)');
xticks([0 0.2 0.4 0.6 0.8 1.0]);
yticks([60 70 80 90 100 110]);
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


% functions
function LHV = LHV_mass(gas)
    % Find the lower heating value of a gas.  Use extra oxygen to ensure
    % complete combustion.  Since To is so low, dissociation is not an issue.
    % Note that Cantera can have trouble with 25C so use 300K.  Ug!
    % C.F. Edwards, 1-10-10

    % Set a fractional tolerance value.
    toler = 1e-6;

    M    = molecularWeights(gas);
    Nsp  = nSpecies(gas);
    Nel  = nElements(gas);
    iO2  = speciesIndex(gas,'O2');
    iH2O = speciesIndex(gas,'H2O');
    iCO2  = speciesIndex(gas,'CO2');
    iSO2  = speciesIndex(gas,'SO2');

    % Save the state so we can reset it at the end.
    Tin = temperature(gas);
    Pin = pressure(gas);
    xin = moleFractions(gas);

    % Set the temperature and pressure at which to evaluate the heating value.
    % Sometimes Cantera will crash Matlab if you use 25C.  So you might want to
    % use 300K instead.  It makes no difference to the heating value.
    % Tref = 25+273.15;   % K
    Tref = 300;         % K
    Pref = 101325;      % Pa
    % Use one mole of fuel mixture.
    Nfuel = moleFractions(gas); 

    % Cantera can have trouble equilibrating with small, nonzero values for
    % trace species.  It will make no difference to the heating value if we
    % kill them.
    for j=1:1:Nsp
        if(Nfuel(j) < toler)
            Nfuel(j) = 0;
        end
    end
    % Renormalize and get the molecular mass.
    Nfuel = Nfuel/sum(Nfuel);
    set(gas,'T',Tref,'P',Pref,'X',Nfuel);
    Mfuel = meanMolecularWeight(gas);

    % Find out which elements the gas contains.  We only care about C,H,O,S.
    Has_H = 0;
    Has_C = 0;
    Has_S = 0;
    Has_O = 0;
    for i=1:1:Nel
        if(strcmp(elementName(gas,i),'H'))
            Has_H = 1;
        elseif(strcmp(elementName(gas,i),'C'))
            Has_C = 1;
        elseif(strcmp(elementName(gas,i),'S'))
            Has_S = 1;
        elseif(strcmp(elementName(gas,i),'O'))
            Has_O = 1;
        end
    end

    % Find the stoichiometric oxygen requirment for this fuel.
    N_H = 0;
    N_C = 0;
    N_S = 0;
    N_O = 0;
    for i=1:1:Nsp
        This_Species_Counts = (Nfuel(i)~=0.0)&&(i~=iH2O)&&(i~=iCO2)&&(i~=iSO2);
        if This_Species_Counts
            % function n = nAtoms(gas,i,m)
            % NATOMS-Number of atoms of m in species k.
            if(Has_H)
                m = elementIndex(gas,'H');
                N_H  = N_H + Nfuel(i)*nAtoms(gas,i,m);
            end
            if(Has_C)
                m = elementIndex(gas,'C');
                N_C  = N_C + Nfuel(i)*nAtoms(gas,i,m);
            end
            if(Has_S)
                m = elementIndex(gas,'S');
                N_S  = N_S + Nfuel(i)*nAtoms(gas,i,m);
            end
            if(Has_O)
                m = elementIndex(gas,'O');
                N_O  = N_O + Nfuel(i)*nAtoms(gas,i,m);
            end
        end
    end
    % We have the total number of each of the relevant atoms in the fuel.  Now
    % find the oxygen needed by the H, C, and S atoms.  Note that if the fuel
    % already contains some oxygen that will be present above and beyond the
    % amount we choose to add.
    O_needed   = N_H/2 + 2*N_C + 2*N_S;
    O2_needed  = O_needed/2;

    % Check to see if this is not combustible.
    if O_needed == 0
        % Restore state of gas.
        set(gas,'T',Tin,'P',Pin,'X',xin);
        % Set the LHV to zero.
        LHV = 0;
        return
    end

    % Add enough oxygen to get 100% excess air.  This will be enough to
    % suppress CO and H2.
    Noxid = 2*O2_needed;
    % Make a mixture.
    Nmix = Nfuel;                   % Put fuel in the mixture
    Nmix(iO2) = Nmix(iO2) + Noxid;  % Add oxygen in excess of existing

    % Find the mass of the mix (fuel plus added oxid).
    mass_mix  = 0;
    for j=1:1:Nsp
        mass_mix  = mass_mix  + Nmix(j)*M(j);
    end
    % Since there is one mole of fuel in the mix, the molecular mass gives
    % its contribution to the total mixture mass.
    mass_fraction_fuel  = Mfuel/mass_mix;
    % Normalize the mix.
    Xmix = Nmix/sum(Nmix);
    % Do the heating value problem.
    set(gas,'T',Tref,'P',Pref,'X',Xmix);        % Put in the mixture
    h_reactants = enthalpy_mass(gas);           % Find its enthalpy
    equilibrate(gas,'TP');                      % Burn it at constant TP
    h_products = enthalpy_mass(gas);            % Find the enthalpy

    % The enthalpy difference per unit mass of fuel is the heating value.
    % Since the water is all vapor, this is the lower heating value.
    LHV = (h_reactants - h_products)/mass_fraction_fuel;

    % Restore state of gas.
    set(gas,'T',Tin,'P',Pin,'X',xin);
end



function x = flowExergy_mass(fluid)
    % Return the specific flow exergy of the fluid (J/kg);

    global To Po mu_o
    global g_tm
        

    % If g_tm did not get declared outside of this routine, force its
    % calculation.
    if isempty(g_tm)
        disp('Global variable g_tm does not exist...')
        g_tm = 0;
    end

    % Get the basics from the fluid.
    h = enthalpy_mass(fluid);
    s = entropy_mass(fluid);
    b = h - To*s;

    % Figure out which fluid you have and deal with it accordingly.
    % Use the first/only species name to check.
    species1 = char(speciesName(fluid,1));

    % Check two character names.
    if((nSpecies(fluid) == 1) && (length(species1) == 2))
        if(strcmp(species1,'H2'))
            if(g_tm == 0)   % Don't calc g_tm if already set
                ref = importPhase('liquidvapor.xml','hydrogen');
                set(ref,'T',To,'P',Po);
                g_tm = gibbs_mass(ref);
            end
            xchem = 1.165e+008;
            x = b - g_tm + xchem;
            return
        end

        if(strcmp(species1,'N2'))
            if(g_tm == 0)   % Don't calc g_tm if already set
                ref = importPhase('liquidvapor.xml','nitrogen');
                set(ref,'T',To,'P',Po);
                g_tm = gibbs_mass(ref);
            end
            xchem = 2.503e+004;
            x = b - g_tm + xchem;
            return
        end

        if(strcmp(species1,'O2'))
            if(g_tm == 0)   % Don't calc g_tm if it is already set
                ref = importPhase('liquidvapor.xml','oxygen');
                set(ref,'T',To,'P',Po);
                g_tm = gibbs_mass(ref);
            end
            xchem = 1.239e+005;
            x = b - g_tm + xchem;
            return
        end
    end

    % Check three character names.
    if((nSpecies(fluid) == 1) && (length(species1) == 3))
        if(strcmp(species1,'H2O'))
            if(g_tm == 0)   % Don't calc g_tm if already set
                ref = importPhase('liquidvapor.xml','water');
                set(ref,'T',To,'P',Po);
                g_tm = gibbs_mass(ref);
            end
            % Assume the dead state is saturated so xchem = 0.
            x = b - g_tm;
            return
        end

        if(strcmp(species1,'CO2'))
            if(g_tm == 0)   % Don't calc g_tm if it is already set
                ref = importPhase('liquidvapor.xml','carbondioxide');
                set(ref,'T',To,'P',Po);
                g_tm = gibbs_mass(ref);
            end
            xchem = 4.456e+005;
            x = b - g_tm + xchem;
            return
        end
    end

    % Check six character names.
    if((nSpecies(fluid) == 1) && (length(species1) == 6))
        if(strcmp(species1,'C2F4H2'))
            if(g_tm == 0)   % Don't calc g_tm if it is already set
                ref = importPhase('liquidvapor.xml','hfc134a');
                set(ref,'T',To,'P',Po);
                g_tm = gibbs_mass(ref);
            end
            % R134a does not react to environmental species, so set xchem = 0.
            x = b - g_tm;
            return
        end
        % If we ever have other 6 char names, put them here!
    end

    % If you get down to here, it is the GRI30 property set so do full chemical
    % exergy.  Note that the dead state chemical potentials must already be
    % loaded into the global array mu_o.  (If you calculated them here,
    % repeated calls would take forever.)

    % When a non-environmental species is converted it does so by the reaction:
    % CxHyOzNwARq -> xCO2 + (y/2)H2O + (w/2)N2 + qAR - (x + y/4 - z/2)O2
    % Note that this also treats environmental species present in the resource
    % correctly--they are converted to themselves by the reaction!

    N  = moleFractions(fluid);
    M  = meanMolecularWeight(fluid);
    Ns = nSpecies(fluid);

    % Define indices to find the environmental species in the array.
    % Note that argon as a species goes by "AR" while argon as an element goes
    % by "Ar".  Isn't Cantera wonderful some times?
    iCO2 = speciesIndex(fluid,'CO2');
    iH2O = speciesIndex(fluid,'H2O');
    iN2  = speciesIndex(fluid,'N2');
    iO2  = speciesIndex(fluid,'O2');
    iAR  = speciesIndex(fluid,'AR');

    G_o = 0;
    for i=1:1:Ns
        if(N(i) ~= 0.0)
            % function n = nAtoms(a,k,m)
            % NATOMS-Number of atoms of m in species k.
            N_O  = nAtoms(fluid,i,elementIndex(fluid,'O')); 
            N_H  = nAtoms(fluid,i,elementIndex(fluid,'H'));
            N_C  = nAtoms(fluid,i,elementIndex(fluid,'C'));
            N_N  = nAtoms(fluid,i,elementIndex(fluid,'N'));
            N_AR = nAtoms(fluid,i,elementIndex(fluid,'Ar'));

            N_CO2 = N_C;
            N_H2O = N_H/2;
            N_N2  = N_N/2;
            N_O2  = N_O/2 - N_C - N_H/4;

            G_o = G_o + N(i)*(N_CO2*mu_o(iCO2) + (N_H2O)*mu_o(iH2O) + (N_N2)*mu_o(iN2) + (N_AR)*mu_o(iAR) + (N_O2)*mu_o(iO2));
        end
    end
    x = b - G_o/M;
end

function x = exergy_mass(fluid)
    % Return the mass-specific internal exergy of the fluid (J/kg);
    % C.F. Edwards, 1-10-10

    % global To Po mu_o
    global g_tm
    To = 298.15;
    Po = 101325;
    mu_o = chemPotentials(fluid);

    % The variable g_tm is used to save the Gibbs function of a pure fluid in
    % the thermomechanical dead state.  This saves time recalculating it.  Note
    % that it will not be correct if you change fluids.  So if you do, reset
    % its value to zero in the main routine to force its recalculation.

    % If g_tm did not get initialized outside of this routine, do so by  
    % forcing its calculation.
    if isempty(g_tm)
        g_tm = 0;
    end

    % Get the basics from the fluid.
    u = intEnergy_mass(fluid);
    r = density(fluid);
    s = entropy_mass(fluid);
    a = u + Po/r - To*s;

    % Figure out which fluid you have and deal with it accordingly.
    % Use the first/only species name to check.
    species1 = char(speciesName(fluid,1));

    % Check two character names.
    if((nSpecies(fluid) == 1) && (length(species1) == 2))
        if(strcmp(species1,'H2'))
            if(g_tm == 0)   % Don't calc g_tm if it is already set
                ref = Solution('liquidvapor.xml','hydrogen');
                set(ref,'T',To,'P',Po);
                g_tm = gibbs_mass(ref);
            end
            xchem = 1.165e+008;
            x = a - g_tm + xchem;
            return
        end

        if(strcmp(species1,'N2'))
            if(g_tm == 0)   % Don't calc g_tm if it is already set
                ref = Solution('liquidvapor.xml','nitrogen');
                set(ref,'T',To,'P',Po);
                g_tm = gibbs_mass(ref);
            end
            xchem = 2.503e+004;
            x = a - g_tm + xchem;
            return
        end

        if(strcmp(species1,'O2'))
            if(g_tm == 0)   % Don't calc g_tm if it is already set
                ref = Solution('liquidvapor.xml','oxygen');
                set(ref,'T',To,'P',Po);
                g_tm = gibbs_mass(ref);
            end
            xchem = 1.239e+005;
            x = a - g_tm + xchem;
            return
        end
    end

    % Check three character names.
    if((nSpecies(fluid) == 1) && (length(species1) == 3))
        if(strcmp(species1,'H2O'))
            if(g_tm == 0)   % Don't calc g_tm if it is already set
                ref = Solution('liquidvapor.xml','water');
                set(ref,'T',To,'P',Po);
                g_tm = gibbs_mass(ref);
            end
            % Assume the dead state is saturated so xchem = 0.
            x = a - g_tm;
            return
        end

        if(strcmp(species1,'CO2'))
            if(g_tm == 0)   % Don't calc g_tm if it is already set
                ref = Solution('liquidvapor.xml','carbondioxide');
                set(ref,'T',To,'P',Po);
                g_tm = gibbs_mass(ref);
            end
            xchem = 4.456e+005;
            x = a - g_tm + xchem;
            return
        end
    end

    % Check six character names.
    if((nSpecies(fluid) == 1) && (length(species1) == 6))
        if(strcmp(species1,'C2F4H2'))
            if(g_tm == 0)   % Don't calc g_tm if it is already set
                ref = Solution('liquidvapor.xml','hfc134a');
                set(ref,'T',To,'P',Po);
                g_tm = gibbs_mass(ref);
            end
            % R134a does not react to environmental species, so set xchem = 0.
            x = a - g_tm;
            return
        end
        % If we ever have other 6 char names, put them here!
    end

    % If you get down to here, it is the GRI30 property set so do full chemical
    % exergy.  Note that the dead state chemical potentials must already be
    % loaded into the global array mu_o.  (If you calculated them here,
    % repeated calls would take forever.)

    % When a non-environmental species is converted it does so by the reaction:
    % CxHyOzNwARq -> xCO2 + (y/2)H2O + (w/2)N2 + qAR - (x + y/4 - z/2)O2
    % Note that this also treats environmental species present in the resource
    % correctly--they are converted to themselves by the reaction!

    N  = moleFractions(fluid);
    M  = meanMolecularWeight(fluid);
    Ns = nSpecies(fluid);

    % Define indices to find the environmental species in the array.
    % Note that argon as a species goes by "AR" while argon as an element goes
    % by "Ar".  Isn't Cantera wonderful some times?
    iCO2 = speciesIndex(fluid,'CO2');
    iH2O = speciesIndex(fluid,'H2O');
    iN2  = speciesIndex(fluid,'N2');
    iO2  = speciesIndex(fluid,'O2');
    iAR  = speciesIndex(fluid,'AR');

    G_o = 0;
    for i=1:1:Ns
        if(N(i) ~= 0.0)
            % function n = nAtoms(a,k,m)
            % NATOMS-Number of atoms of m in species k.
            N_O  = nAtoms(fluid,i,elementIndex(fluid,'O'));
            N_H  = nAtoms(fluid,i,elementIndex(fluid,'H'));
            N_C  = nAtoms(fluid,i,elementIndex(fluid,'C'));
            N_N  = nAtoms(fluid,i,elementIndex(fluid,'N'));
            N_AR = nAtoms(fluid,i,elementIndex(fluid,'Ar'));

            N_CO2 = N_C;
            N_H2O = N_H/2;
            N_N2  = N_N/2;
            N_O2  = N_O/2 - N_C - N_H/4;

            G_o = G_o + N(i)*(N_CO2*mu_o(iCO2) + (N_H2O)*mu_o(iH2O) + (N_N2)*mu_o(iN2)...
                + (N_AR)*mu_o(iAR) + (N_O2)*mu_o(iO2));
        end
    end
    x = a - G_o/M;
end
