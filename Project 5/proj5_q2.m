% Gasification using ideal gases in gasification_small.xml
% Does the problem backwards and uses general HCs or oxygenates.
% C.F. Edwards, 1/25/09

clear all
format compact
fprintf('\n************************************************************\n')

% Input Conditions
Tin = 150 + 273.15;
Pin = 55*oneatm;

% Get the species.
gas  = Solution('gasification_small.xml');
Nsp  = nSpecies(gas);
Ma   = atomicMasses(gas);
Mm   = molecularWeights(gas);
iCH4 = speciesIndex(gas,'CH4');
iO2  = speciesIndex(gas,'O2');
iN2  = speciesIndex(gas,'N2');
iH2O = speciesIndex(gas,'H2O');
iCO  = speciesIndex(gas,'CO');
iCO2 = speciesIndex(gas,'CO2');
iH2  = speciesIndex(gas,'H2');
iSO2 = speciesIndex(gas,'SO2');
iC   = elementIndex(gas,'C');
iH   = elementIndex(gas,'H');
iO   = elementIndex(gas,'O');
iN   = elementIndex(gas,'N');
iS   = elementIndex(gas,'S');
iAr  = elementIndex(gas,'AR');

% Gasify methane using air.  Use molar ratios to specify the oxygen and water to be
% used.  For complete combustion CH4 + 2O2 -> CO2 + 2H2O so we know that we
% want an O2/C ratio of less than 2:1.  We will add water to keep the
% temperature down and to provide some WGS action.  Normally the values are
% around H2O/C of 1:1.

% Operation matrix
H2OC = 1.5:-0.01:0;
n_h2o = size(H2OC,2);
O2C=0.5:-0.01:0.1;
n_o2c = size(O2C,2);

% Initialize Data Storage
Temp_data(n_o2c,n_h2o) = 0;
Species_data(n_o2c,n_h2o,Nsp) = 0;
CGE_data(n_o2c,n_h2o) = 0;
Syngas_yield_data(n_o2c,n_h2o) = 0;


% Now do the problem backwards.  
% Make a loop to vary the oxygen-carbon ratio.

for j = 1:n_h2o
    H2O_C = H2OC(j); % water to carbon ratio 
for i = 1:n_o2c 
    O2_C = O2C(i);  % oxygen to carbon ratio
    
    % Number of atoms from Ultimate analysis
    Nin1_coal   = 1;
    num_C_coal  = Nin1_coal * 0.755/Ma(iC);
    num_H_coal  = Nin1_coal * 0.064/Ma(iH);
    num_O_coal  = Nin1_coal * 0.152/Ma(iO);
    num_N_coal  = Nin1_coal * 0.015/Ma(iN);
    num_S_coal  = Nin1_coal * 0.014/Ma(iS);
    
    % number of atoms of 1 mol of C in Coal
    num_C     = num_C_coal/num_C_coal;
    num_H     = num_H_coal/num_C_coal;
    num_O     = num_O_coal/num_C_coal;
    num_N     = num_N_coal/num_C_coal;
    num_S     = num_S_coal/num_C_coal;

    % Proximate analysis
    Fix_C     = 0.30;   
    Vol_matter = 0.23;
    Water     = 0.07;
    Ash       = 0.40;
    
    % Inlet H2O and Air
    Nin1_O2   = O2_C * num_C;
    Nin1_H2O  = H2O_C * num_C;
    Nin1_N2   = 3.76 * Nin1_O2;
    N2_C      = 3.76 * O2_C;

    M_maf     = num_C*Ma(iC) + num_H*Ma(iH) + num_S*Ma(iS) + num_O*Ma(iO) + num_N*Ma(iN);  % [kg/kmol]
    
    % Mass Flow Rate of Coal
    mdot_maf = 1*M_maf; % [kmol/s] * [kg/kmol] = [kg/s] of coal
    mdot_ash = mdot_maf*Ash/(Fix_C+Vol_matter);  % kg/s of ash per kg of maf coal
    mdot_mf  = mdot_maf*(Fix_C+Vol_matter+Ash)/(Fix_C+Vol_matter); % kg/s of mf coal
    mdot_slag= mdot_ash; % assuming all ash turns into slag
    ash_maf  = mdot_slag/mdot_maf; % kg of slag per kg of maf
    mdot_mois= mdot_maf*Water/(Fix_C+Vol_matter); % mol of water in af coal per second
    ndot_mois= mdot_mois/18; % moles of water per second

    % IL #6 Ash content:  Sio2: 47.5, Al2O3: 17.9, TiO2: 0.8, Fe2O3: 20.1, 
    % CaO:5.8, MgO: 1, Na2O: 0.4, K2O: 1.8, SO3: 4.6, P2O5: 0.1
    ash_content = [47.5, 17.9, 0.8, 20.1, 5.8, 1, 0.4, 1.8, 4.6, 0.1]*1e-2;
    ash_species_mm = [60 102 80 160 56 40.3 62 94.2 80.6 283.9];
    ash_o_species = [2 3 2 3 1 1 1 1 3 5];
    mdot_o_ash =  sum(ash_o_species*16.*ash_content) / sum(ash_content.*ash_species_mm); % kg oxygen per kg of ash for slag to oxidize
    ndot_o2_ash = mdot_o_ash/32; % mol of O2 per kg of ash

    % number of mole of product gases
    Nin2        = zeros(Nsp,1);
    Nin2(iCO)   = 1 * num_C;
    Nin2(iH2)   = 0.5 * num_H;
    Nin2(iH2O)  = Nin1_H2O;
    Nin2(iSO2)  = 1 * num_S;
    Nin2(iN2)   = Nin1_N2 + 1*0.5*num_N;
    Nin2(iO2)   = Nin1_O2 + 1*0.5*num_O - 0.5*Nin2(iCO) - 1*Nin2(iSO2);
    if(Nin2(iO2) < 0)
        Nin2(iH2O) = Nin2(iH2O) - (-1)*Nin2(iO2)*2;
        Nin2(iH2)  = Nin2(iH2)  + (-1)*Nin2(iO2)*2;
        Nin2(iO2)  = 0;
    end
   

    % Use Cantera to find the enthalpies of the product gases for the
    % heating value problem.
    set(gas,'T',25+273.15,'P',oneatm,'X','CO2:1');
    hbar_CO2 = enthalpy_mole(gas);
    set(gas,'T',25+273.15,'P',oneatm,'X','H2O:1');
    hbar_H2O = enthalpy_mole(gas);
    set(gas,'T',25+273.15,'P',oneatm,'X','SO2:1');
    hbar_SO2 = enthalpy_mole(gas);
    
    % TODO: contest this LHV inversion

    % Data for input Coal             
    Qlhv = 32.1e6;                 % LHV is known for input coal, J/kg
    Cp    = 0.7e3;                  % Cp is known, J/kg-K

    % Finish Inverting LHV Problem: Enthalpy for input coal
    h_maf = num_C*hbar_CO2/M_maf + num_H/2*hbar_H2O/M_maf + num_S*hbar_SO2/M_maf + Qlhv;

    % Add in energy from heating coal to Tin, use the average Cp. 
    h_mf_delta = Cp*(Tin-(225+273.1598));
      
    % Find specific enthalpies of air and water inputs
    set(gas,'T',Tin,'P',oneatm,'X','H2O:1');
    h_H2O_Tin = enthalpy_mole(gas);
    set(gas,'T',Tin,'P',oneatm,'X','O2:1');
    h_O2_Tin = enthalpy_mole(gas);
    set(gas,'T',Tin,'P',oneatm,'X','N2:1');
    h_N2_Tin = enthalpy_mole(gas);
    
    %% Ash and Slag Properties
    %Tm      = 1350; % [oC]
    Cp_slag = 1e3; % J/kg*K
    Cp_ash  = 1e3; % J/kg*K
    hfg_ash = 240e3; % J/kg
    reducing_T = 1270;  % [oC] % TODO: understand what these are
    oxidizing_T = 1430; % [oC]
    
    %% Find Equillibrium temperature Iteratively

    T_predicted_min=400; % guess
    T_predicted_max=2500;
    while (1)
        T_predicted_avg = (T_predicted_max + T_predicted_min)/2; % bisection search
        
        % Reset H2 and H2O values
        Nin2(iH2)   = 0.5 * num_H;
        Nin2(iH2O)  = Nin1_H2O;

        % Account for ash/slag
        if T_predicted_avg<reducing_T    % ash not melted
            Nin2(iO2)  = Nin1_O2 + 1*0.5*num_O - 0.5*Nin2(iCO) - 1*Nin2(iSO2);
            h_slag_delta = Cp_ash*(T_predicted_avg-25); 
        elseif T_predicted_avg>reducing_T && T_predicted_avg<oxidizing_T % ash melted
            Nin2(iO2)   = Nin1_O2 + 1*0.5*num_O - 0.5*Nin2(iCO) - 1*Nin2(iSO2);
            h_slag_delta = Cp_ash*(reducing_T-25) + hfg_ash + Cp_slag*(T_predicted_avg-reducing_T); 
        else % ash melted and oxidized
            h_slag_delta = Cp_ash*(reducing_T-25) + hfg_ash + Cp_slag*(T_predicted_avg-reducing_T); 
            % Adjust for oxidized slag
            Nin2(iO2)   = Nin1_O2 + 1*0.5*num_O - 0.5*Nin2(iCO) - 1*Nin2(iSO2) - ndot_o2_ash*ash_maf;
        end

        % Check oxygen value isn't negative
        if(Nin2(iO2) < 0)
                Nin2(iH2O) = Nin2(iH2O) - (-1)*Nin2(iO2)*2;
                Nin2(iH2)  = Nin2(iH2)  + (-1)*Nin2(iO2)*2;
                Nin2(iO2)  = 0;
        end

        mdot_syngas   = Nin2(iCO)*Mm(iCO) + Nin2(iH2)*Mm(iH2) + Nin2(iH2O)*Mm(iH2O) + Nin2(iSO2)*Mm(iSO2) + Nin2(iN2)*Mm(iN2) + Nin2(iO2)*Mm(iO2);

        h_syngas = mdot_maf/mdot_syngas * (h_maf + (1+ash_maf)*h_mf_delta + O2_C*h_O2_Tin/M_maf + N2_C*O2_C*h_N2_Tin/M_maf + (H2O_C)*h_H2O_Tin/M_maf - (ash_maf)*h_slag_delta); %%%%
    
        % Use the syngas composition and backwards-computed enthalpy to find the
        % equilibrium state.
        set(gas,'T',300,'P',oneatm,'X',Nin2);   % Pre-equilibrate at low T
        equilibrate(gas,'TP');                  % This can soothe a cranky Cantera.
        set(gas,'H',h_syngas,'P',Pin);  % Set the correct state.
        equilibrate(gas,'HP');                  % Equilibrate for real.
        ToutC  = temperature(gas) - 273.15;
        xout = moleFractions(gas);
        

        if abs(ToutC-T_predicted_avg)>0.1 % Not within error
            if ToutC>T_predicted_avg
                T_predicted_min = T_predicted_avg; % Change max and min values accordingly
            else
                T_predicted_max = T_predicted_avg;
            end
        else
            Temp_data(i,j) = ToutC;
            Species_data(i,j,:) = xout;
            %CGE_data TODO: write these metrics
            %Syngas_yield_data
            break % T_predicted_avg within 0.1 C, exit loop
        end
    end

end
   fprintf("Iteration %d out of %d \n", j, n_h2o)   
end

%% PLOTTING CODE
