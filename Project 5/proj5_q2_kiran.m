% Project 5 Part 2

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
iH2S = speciesIndex(gas,'H2S');
iNH3 = speciesIndex(gas,'NH3');
iCOS = speciesIndex(gas,'COS');
iC   = elementIndex(gas,'C');
iH   = elementIndex(gas,'H');
iO   = elementIndex(gas,'O');
iN   = elementIndex(gas,'N');
iS   = elementIndex(gas,'S');

% Gasify methane using air.  Use molar ratios to specify the oxygen and water to be
% used.  For complete combustion CH4 + 2O2 -> CO2 + 2H2O so we know that we
% want an O2/C ratio of less than 2:1.  We will add water to keep the
% temperature down and to provide some WGS action.  Normally the values are
% around H2O/C of 1:1.

% Operation matrix
H2OC = 1.5:-0.02:0;
n_h2o = size(H2OC,2);
O2C=0.5:-0.02:0.1;
n_o2c = size(O2C,2);

% Initialize Data Storage
Temp_data(n_o2c,n_h2o) = 0;
Species_data(n_o2c,n_h2o,Nsp) = 0;
% CGE_data(n_o2c,n_h2o) = 0;
% Syngas_yield_data(n_o2c,n_h2o) = 0;
ratio_C(n_o2c,n_h2o) = 0;
ratio_coal(n_o2c,n_h2o) = 0;
mean_weight(n_o2c,n_h2o) = 0;
CGE(n_o2c,n_h2o) = 0;
Syngas_yield(n_o2c,n_h2o) = 0;

xtest = zeros(16,1);

% Now do the problem backwards.  
% Make a loop to vary the oxygen-carbon ratio.

for j = 1:n_h2o
    H2O_C = H2OC(j); % water to carbon ratio 
for i = 1:n_o2c 
    O2_C = O2C(i);  % oxygen to carbon ratio
    
    % Number of atoms from Ultimate analysis
    Nin1_coal   = 1;
    num_C_coal  = Nin1_coal * 0.784/Ma(iC);
    num_H_coal  = Nin1_coal * 0.054/Ma(iH);
    num_O_coal  = Nin1_coal * 0.099/Ma(iO);
    num_N_coal  = Nin1_coal * 0.014/Ma(iN);
    num_S_coal  = Nin1_coal * 0.049/Ma(iS);
    
    % number of atoms per 1 mol of C in Coal
    num_C     = num_C_coal/num_C_coal;
    num_H     = num_H_coal/num_C_coal;
    num_O     = num_O_coal/num_C_coal;
    num_N     = num_N_coal/num_C_coal;
    num_S     = num_S_coal/num_C_coal;

    % Proximate analysis
    Fix_C     = 0.393;   
    Vol_matter = 0.37;
    Water     = 0.13;
    Ash       = 0.107;
    
    % Inlet H2O and Air
    Nin1_O2   = O2_C * num_C;
    Nin1_H2O  = H2O_C * num_C;
    Nin1_N2   = 0; % 3.76 * Nin1_O2; using neat oxygen
    N2_C      = 0; % 3.76 * O2_C; using neat oxygen

    M_maf     = num_C*Ma(iC) + num_H*Ma(iH) + num_S*Ma(iS) + num_O*Ma(iO) + num_N*Ma(iN);  % [kg maf/ kmol C]
    
    % Mass Flow Rate of Coal
    mdot_maf = 1*M_maf; % [kmol C/s] * [kg maf/kmol C] = [kg maf/s] of coal
    mdot_ash = mdot_maf*Ash/(Fix_C+Vol_matter);  % kg/s of ash per kg of maf coal
    mdot_mf  = mdot_maf*(Fix_C+Vol_matter+Ash)/(Fix_C+Vol_matter); % kg/s of mf coal
    mdot_slag= mdot_ash; % assuming all ash turns into slag
    ash_maf  = mdot_slag/mdot_maf; % kg of slag per kg of maf
    mdot_mois= mdot_maf*Water/(Fix_C+Vol_matter); % mol of water in af coal per second
    ndot_mois= mdot_mois/18; % moles of water per second

    % number of mole of product gases
    Nin2        = zeros(Nsp,1);
    Nin2(iCO)   = 1 * num_C; 
    Nin2(iH2)   = 0.5*num_H + Nin1_H2O; % mol H2/ mol C
    Nin2(iSO2)  = 1 * num_S;
    Nin2(iN2)   = Nin1_N2 + 0.5*num_N;
    Nin2(iO2)   = Nin1_O2 + 0.5*num_O + 0.5*(Nin1_H2O) - 0.5*Nin2(iCO) - 1*Nin2(iSO2);
   
    flag = 0; % Pyrolysis Condition
    if(Nin2(iO2) < 0)
        flag = 1;
        Nin2(iO2) = 0;
    end
   
    % Unit [kg/mol C]
    mdot_syngas   = Nin2(iCO)*Mm(iCO) + Nin2(iH2)*Mm(iH2) + Nin2(iSO2)*Mm(iSO2) + Nin2(iN2)*Mm(iN2) + Nin2(iO2)*Mm(iO2);

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
    Qlhv = 33.7e6;                 % LHV is known for input coal, J/kg
    Cp    = 0.7e3;                  % Cp is known, J/kg-K

    % Finish Inverting LHV Problem: Enthalpy for input coal
    h_maf = num_C*hbar_CO2/M_maf + num_H/2*hbar_H2O/M_maf + num_S*hbar_SO2/M_maf + Qlhv;

    % Add in energy from heating coal to Tin, use the average Cp. 
    h_mf_delta = Cp*(Tin-(25+273.1598));
      
    % Find specific enthalpies of air and water inputs
    set(gas,'T',Tin,'P',oneatm,'X','H2O:1');
    h_H2O_Tin = enthalpy_mole(gas);
    set(gas,'T',Tin,'P',oneatm,'X','O2:1');
    h_O2_Tin = enthalpy_mole(gas);
    set(gas,'T',Tin,'P',oneatm,'X','N2:1');
    h_N2_Tin = enthalpy_mole(gas);
    
    % Find Equillibrium temperature Iteratively

    T_predicted_min=100; % guess
    T_predicted_max=3500;

    Nin2 = Nin2./sum(Nin2); % Normalize 

    while (1)
        T_predicted_avg = (T_predicted_max + T_predicted_min)/2; % bisection search

        %h_syngas = (mdot_maf/mdot_syngas)*(h_maf + (1+ash_maf)*h_mf_delta + Nin1_O2*h_O2_Tin/M_maf + Nin1_N2*h_N2_Tin/M_maf + (Nin1_H2O)*h_H2O_Tin/M_maf ); %%%%
        h_syngas = (mdot_maf/mdot_syngas) * (h_maf + (1+ash_maf)*h_mf_delta + O2_C*h_O2_Tin/M_maf + N2_C*h_N2_Tin/M_maf + (H2O_C)*h_H2O_Tin/M_maf); %%%%@@@@@@@@@@@@22
    
        if flag % pyrolysis condition
            ratio_C(i,j) = 0;
            ratio_coal(i,j) = 0;
            mean_weight(i,j) = 0;
            Temp_data(i,j) = 0;
            Species_data(i,j,:) = 0;
            break 
        end

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
            ratio_C(i,j) = Nin2(iCO); % mole fraction of carbon in gas
            
            ratio_coal(i,j) = (mdot_maf/mdot_syngas);
            mean_weight(i,j) = meanMolecularWeight(gas);
            Temp_data(i,j) = ToutC;
            Species_data(i,j,:) = xout;
            xtest = xout;
            break % T_predicted_avg within 0.1 C, exit loop
        end
    end
    
    if flag
        CGE(i,j) = 0;
        Syngas_yield(i,j) = 0;
    else
        CGE(i,j) = mdot_syngas*LHV_mass(xtest)/(mdot_maf*Qlhv)*100;
        MW_Syngas = Mm'*xtest;
        Syngas_yield(i,j) = mdot_syngas/MW_Syngas*(xtest(iH2)+xtest(iCO))/1;
    end
    

end
   fprintf("Iteration %d out of %d \n", j, n_h2o)   
end


% Generate CGE and SYngas Yield
set(gas,'T',Tin,'P',oneatm,'X','N2:1');
h_N2_Tin = enthalpy_mole(gas);

% for i = 1:n_o2c
%     for j = 1:n_h2o 
%         co = Species_data(i,j,iCO);
%         h2 = Species_data(i,j,iH2);
%         % Syngas_yield_data(i,j) = (h2 + co)/ratio_C(i,j);
%         % heating values: 120 MJ/kg for hydrogen, 10.16 MJ/kg for CO
%         %CGE_data(i,j) = 1e2*(120e6*h2*2 + 10.16e6*co*28) /(Qlhv*ratio_coal(i,j)*mean_weight(i,j));% Fix with enthalpies fo hydrogen and CO
%     end
% end

% PLOTTING CODE\
figure(1)
contour(O2C, H2OC, Temp_data', 600:100:2500,'ShowText','on', 'LineWidth',2)
xlabel("Oxygen-Carbon Ratio")
ylabel("Water-Carbon Ratio")
title("Temperature (C)")
%improvePlot

figure(2)
contour(O2C, H2OC, Species_data(:,:,iH2)','ShowText','on', 'LineWidth',2)
xlabel("Oxygen-Carbon Ratio")
ylabel("Water-Carbon Ratio")
title("H2 Mole Ratio")
%improvePlot

figure(3)
contour(O2C, H2OC, CGE(:,:)','ShowText','on', 'LineWidth',2)
xlabel("Oxygen-Carbon Ratio")
ylabel("Water-Carbon Ratio")
title("Cold Gas Efficiency (%)")
%improvePlot

figure(4)
contour(O2C, H2OC, Syngas_yield(:,:)','ShowText','on', 'LineWidth',2)
xlabel("Oxygen-Carbon Ratio")
ylabel("Water-Carbon Ratio")
title("Syngas Molar Yield")
%improvePlot

figure(5)
contour(O2C, H2OC, Species_data(:,:,iCO)','ShowText','on', 'LineWidth',2)
xlabel("Oxygen-Carbon Ratio")
ylabel("Water-Carbon Ratio")
title("CO Mole Ratio")

figure(6)
contour(O2C, H2OC, Species_data(:,:,iCH4)','ShowText','on', 'LineWidth',2)
xlabel("Oxygen-Carbon Ratio")
ylabel("Water-Carbon Ratio")
title("CH_4 Mole Ratio")

figure(7)
contour(O2C, H2OC, Species_data(:,:,iH2O)','ShowText','on', 'LineWidth',2)
xlabel("Oxygen-Carbon Ratio")
ylabel("Water-Carbon Ratio")
title("H_2O Mole Ratio")

figure(8)
contour(O2C, H2OC, Species_data(:,:,iN2)','ShowText','on', 'LineWidth',2)
xlabel("Oxygen-Carbon Ratio")
ylabel("Water-Carbon Ratio")
title("N_2 Mole Ratio")

figure(9)
contour(O2C, H2OC, Species_data(:,:,iCO2)','ShowText','on', 'LineWidth',2)
xlabel("Oxygen-Carbon Ratio")
ylabel("Water-Carbon Ratio")
title("CO_2 Mole Ratio")

figure(10)
contour(O2C, H2OC, Species_data(:,:,iH2S)','ShowText','on', 'LineWidth',2)
xlabel("Oxygen-Carbon Ratio")
ylabel("Water-Carbon Ratio")
title("H_2S Mole Ratio")

figure(11)
contour(O2C, H2OC, Species_data(:,:,iCOS)','ShowText','on', 'LineWidth',2)
xlabel("Oxygen-Carbon Ratio")
ylabel("Water-Carbon Ratio")
title("COS Mole Ratio")

figure(12)
contour(O2C, H2OC, Species_data(:,:,iNH3)','ShowText','on', 'LineWidth',2)
xlabel("Oxygen-Carbon Ratio")
ylabel("Water-Carbon Ratio")
title("NH3 Mole Ratio")

function LHV_mass=LHV_mass(x_gas)
    T=25+273.15; % [K]
    P=101325; % [Pa]     

    gas = Solution('gasification_small.xml');
    nsp = nSpecies(gas);
    set(gas,'P',P,'T',T,'X',x_gas);
    M=molecularWeights(gas); % molecular weight vector
    mass_fuel=x_gas'*M;  % mass of 1 mol of fuel 

    iCH4  = speciesIndex(gas,'CH4'); iH2  = speciesIndex(gas,'H2'); iCO  = speciesIndex(gas,'CO');
    iH2O  = speciesIndex(gas,'H2O'); iCO2  = speciesIndex(gas,'CO2'); iNH3  = speciesIndex(gas,'NH3');
    iO2   = speciesIndex(gas,'O2'); iNO  = speciesIndex(gas,'NO'); iNO2  = speciesIndex(gas,'NO2'); iN2O  = speciesIndex(gas,'N2O');
    iN2   = speciesIndex(gas,'N2'); iAR  = speciesIndex(gas,'AR'); iH2S  = speciesIndex(gas,'H2S');
    iSO2  = speciesIndex(gas,'SO2'); iCOS  = speciesIndex(gas,'COS'); iHCN  = speciesIndex(gas,'HCN');

    % number of atoms
    nO=0; nH=0; nC=0; nN=0; nAR=0; nS=0;
    for i=1:nsp
        nO = nO + nAtoms(gas,i,elementIndex(gas,'O'));
        nH = nH + nAtoms(gas,i,elementIndex(gas,'H'));
        nC = nC + nAtoms(gas,i,elementIndex(gas,'C'));
        nN = nN + nAtoms(gas,i,elementIndex(gas,'N'));
        nAR = nAR + nAtoms(gas,i,elementIndex(gas,'Ar'));
        nS = nS + nAtoms(gas,i,elementIndex(gas,'S'));
    end

    % composition of fuel and excessive air mixture
    nO_for_stoi=0.5*(2*nC+0.5*nH-1*nO);
    if nC==0 && nH==0
        LHV_mass=0;
    else
        Excess_air=3.0; % 3 times excessive air 
        x_mix=x_gas;
        x_mix(iO2)=x_gas(iO2)+(1+Excess_air)*nO_for_stoi;
        x_mix(iN2)=x_gas(iN2)+(1+Excess_air)*3.76*nO_for_stoi;
        %x_mix=x_mix/sum(x_mix);
    
        % find LHV
        set(gas,'P',P,'T',T,'X',x_mix);
        h_react=enthalpy_mass(gas);
        equilibrate(gas,'TP');
        h_prod=enthalpy_mass(gas);
        LHV=h_react-h_prod; % [J/kg_mixture]
        
        % find LHV [J/kg_fuel]
        mass_mix=x_mix'*M;
        M_fraction_fuel=mass_fuel/mass_mix;
        LHV_mass=LHV/M_fraction_fuel;
    end
end