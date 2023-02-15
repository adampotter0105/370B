% Gasification using ideal gases in gasification_small.xml
% Does the problem backwards and uses general HCs or oxygenates.
% C.F. Edwards, 1/25/09

clear all
format compact
fprintf('\n************************************************************\n')

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
iCOS = speciesIndex(gas,'COS');
iNH3 = speciesIndex(gas,'NH3');
iH2S = speciesIndex(gas,'H2S');

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
H2O_C  = 1;      % Water to carbon ratio

% Now do the problem backwards.  
% Make a loop to vary the oxygen-carbon ratio.

i = 1;
for O2_C=1:-0.01:0 % higher to lower
    O2C(i) = O2_C;  % Oxygen to carbon ratio
    
    % number of atoms from Ultimate analysis
    Nin1_coal   = 1;
    num_C_coal  = Nin1_coal * 0.755/Ma(iC);
    num_H_coal  = Nin1_coal * 0.064/Ma(iH);
    num_O_coal  = Nin1_coal * 0.152/Ma(iO);
    num_N_coal  = Nin1_coal * 0.015/Ma(iN);
    num_S_coal  = Nin1_coal * 0.014/Ma(iS);
    
    % number of atoms of 1 mol of C in Natural Gas
    num_C     = num_C_coal/num_C_coal;
    num_H     = num_H_coal/num_C_coal;
    num_O     = num_O_coal/num_C_coal;
    num_N     = num_N_coal/num_C_coal;
    num_S     = num_S_coal/num_C_coal;

    % Proximate analysis
    Fix_C     = 0.30;   
    Vol_mater = 0.23;
    Water     = 0.07;
    Ash       = 0.40;
    
    % Inlet H2O and Air
    Nin1_O2   = O2_C * num_C;
    Nin1_H2O  = H2O_C * num_C;
    Nin1_N2   = 3.76 * Nin1_O2;
    N2_C      = 3.76 * O2_C;

    M_maf     = num_C*Ma(iC) + num_H*Ma(iH) + num_S*Ma(iS) + (num_O)*Ma(iO) + num_N*Ma(iN);  % [kg/kmol]
    
    %M_ash    = mor_C*Ma(iO);  
    mdot_maf = 1*M_maf; % M_maf = 1*mdot_maf [kg/s]
    mdot_ash = mdot_maf*Ash/(Fix_C+Vol_mater); 
    mdot_mf  = mdot_maf*(Fix_C+Vol_mater+Ash)/(Fix_C+Vol_mater);
    mdot_slag= mdot_ash;
    ash_maf  = mdot_slag/mdot_maf;
    mdot_mois= mdot_maf*Water/(Fix_C+Vol_mater);
    ndot_mois= mdot_mois/18;

    % number of mole of product gases
    Nin2        = zeros(Nsp,1);
    Nin2(iCO)   = 1 * num_C;
    Nin2(iH2)   = 0.5 * num_H;
    Nin2(iH2O)  = Nin1_H2O;         %%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ H2O, No moisture from coal, only feed water
    Nin2(iSO2)  = 1 * num_S;
    Nin2(iN2)   = Nin1_N2 + 1*0.5*num_N;
    Nin2(iO2)   = Nin1_O2 + 1*0.5*num_O - 0.5*Nin2(iCO) - 1*Nin2(iSO2);
    if(Nin2(iO2) < 0)
        Nin2(iH2O) = Nin2(iH2O) - (-1)*Nin2(iO2)*2;
        Nin2(iH2)  = Nin2(iH2)  + (-1)*Nin2(iO2)*2;
        Nin2(iO2)  = 0;
    end
   
    mdot_syngas   = Nin2(iCO)*Mm(iCO) + Nin2(iH2)*Mm(iH2) + Nin2(iH2O)*Mm(iH2O) + Nin2(iSO2)*Mm(iSO2) + Nin2(iN2)*Mm(iN2) + Nin2(iO2)*Mm(iO2);  %@@@@@@@@@@@@@@@@
    Min2_syngas   = mdot_syngas;
    
    % Use Cantera to find the enthalpies of the product gases for the
    % heating value problem.
    set(gas,'T',25+273.15,'P',101325,'X','CO2:1');
    hbar_CO2 = enthalpy_mole(gas);
    set(gas,'T',25+273.15,'P',101325,'X','H2O:1');
    hbar_H2O = enthalpy_mole(gas);
    set(gas,'T',25+273.15,'P',101325,'X','SO2:1');
    hbar_SO2 = enthalpy_mole(gas);
    
    % Here are the data that are available for typical HCs.
    % These are from Cengel and Boles (CB).
                      
    Qlhv_maf = 32.1e6;                 % LHV is known, J/kg
    Cp_mf    = 0.7e3;                  % Cp is known, J/kg-K
    Q_loss   = 0.02*Qlhv_maf;          % loss, J/kg 2%
    
    Qlhv  = Qlhv_maf;
    Cp    = Cp_mf;
    h_maf = num_C*hbar_CO2/M_maf + num_H/2*hbar_H2O/M_maf + num_S*hbar_SO2/M_maf + Qlhv;
    % Use the average Cp.  It's usually all you have.
    h_mf_delta = Cp*(450-25);
      
    % Now that we have the %%%%mole%%%%%-specific enthalpy of the fuel, sum this with
    % the air and water to get the actual enthalpy of the mixture at 450C.
    set(gas,'T',450+273.15,'P',101325,'X','H2O:1');
    h_H2O_450 = enthalpy_mole(gas);
    set(gas,'T',450+273.15,'P',101325,'X','O2:1');
    h_O2_450 = enthalpy_mole(gas);
    set(gas,'T',450+273.15,'P',101325,'X','N2:1');
    h_N2_450 = enthalpy_mole(gas);
    
    Tm      = 1350; % [oC]
    Cp_slag = 1e3;
    Cp_ash  = 1e3;
    hfg_ash = 240e3;
    
    T_predicted_min=500; % guess
    T_predicted_max=2500;
    while (1)
        T_predicted_avg = (T_predicted_max + T_predicted_min)/2;
        if T_predicted_avg<Tm         %%%%%%%%%% T
            h_slag_delta = Cp_ash*(T_predicted_avg-25) + hfg_ash;  %%% ???
        else
            h_slag_delta = Cp_ash*(Tm-25) + hfg_ash + Cp_slag*(T_predicted_avg-Tm);
        end
    
        h_syngas         = mdot_maf/mdot_syngas * (h_maf + (1+ash_maf)*h_mf_delta + O2_C*h_O2_450/M_maf + N2_C*h_N2_450/M_maf + (H2O_C)*h_H2O_450/M_maf - (ash_maf)*h_slag_delta); %%%%@@@@@@@@@@@@22
    
        % Use the syngas composition and backwards-computed enthalpy to find the
        % equilibrium state.
        set(gas,'T',300,'P',100000,'X',Nin2);   % Pre-equilibrate at low T
        equilibrate(gas,'TP');                  % This can soothe a cranky Cantera.
        set(gas,'H',h_syngas,'P',10*100000);  % Set the correct state.
        equilibrate(gas,'HP');                  % Equilibrate for real.
        ToutC(i)  = temperature(gas) - 273.15;
        xout(:,i) = moleFractions(gas);
        

        if abs(ToutC(i)-T_predicted_avg)>0.1
            if ToutC(i)>T_predicted_avg
                T_predicted_min = T_predicted_avg;
            else
                T_predicted_max = T_predicted_avg;
            end
        else
            break
        end
    end

    T_predicted_min=500; % guess
    T_predicted_max=2500;
    while (1)
        T_predicted_avg = (T_predicted_max + T_predicted_min)/2;
        if T_predicted_avg<Tm         %%%%%%%%%% T
            h_slag_delta = Cp_ash*(T_predicted_avg-25) + hfg_ash;
        else
            h_slag_delta = Cp_ash*(Tm-25) + hfg_ash + Cp_slag*(T_predicted_avg-Tm);
        end

        h_syngas_loss     = mdot_maf/mdot_syngas * (h_maf + (1+ash_maf)*h_mf_delta + O2_C*h_O2_450/M_maf + N2_C*h_N2_450/M_maf + H2O_C*h_H2O_450/M_maf - (ash_maf)*h_slag_delta - Q_loss); %%%%
        set(gas,'T',300,'P',100000,'X',Nin2);   % Pre-equilibrate at low T
        equilibrate(gas,'TP');                  % This can soothe a cranky Cantera.
        set(gas,'H',h_syngas_loss,'P',10*100000);  % Set the correct state.
        equilibrate(gas,'HP');  
        ToutC_loss(i)  = temperature(gas) - 273.15;
        xout_loss(:,i) = moleFractions(gas);
        if abs(ToutC(i)-T_predicted_avg)>0.1
            if ToutC(i)>T_predicted_avg
                T_predicted_min = T_predicted_avg;
            else
                T_predicted_max = T_predicted_avg;
            end
        else
            break
        end    
    end

    CGE(i)          = mdot_syngas*LHV_mass(xout_loss(:,i))/(mdot_maf*Qlhv_maf)*100;

    % mdot_maf + mdot_O2 + mdot_N2 + mdot_H2O = mdot_Syngas
    % ndot_syngas = mdot_syngas/MW_Syngas
    % ndot_maf = 1
    MW_Syngas       = Mm'*xout_loss(:,i);
    Syngas_yield(i) = mdot_syngas/MW_Syngas*(xout_loss(iH2,i)+xout_loss(iCO,i))/1;


    i
    i = i+1;
end



figure(1)
clf
hold on
plot(O2C,ToutC,'b')
plot(O2C,ToutC_loss,'b--')
% legend('Methane','Methanol','Gasoline')
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Equilibrium Temperature (\circC)')
text(.05,1800,sprintf('Autothermal Reforming'))
text(.05,1650,'Conditions: 10 bar, 450\circC preheat')
text(.05,1500,sprintf('Water/Carbon Mole Ratio:  %.1f',H2O_C))
scale = axis;
axis([0 1 0 2000])
hold off
plotfixer

figure(2)
clf
hold on
% plot(O2C,xout(iH2,:),'b')
% plot(O2C,xout(iH2O,:),'g')
% plot(O2C,xout(iCH4,:),'r')
% plot(O2C,xout(iCO,:),'m')
% plot(O2C,xout(iCO2,:),'c')
% plot(O2C,xout(iO2,:),'y')
plot(O2C,xout_loss(iH2,:),'b--')
plot(O2C,xout_loss(iH2O,:),'g--')
plot(O2C,xout_loss(iCH4,:),'r--')
plot(O2C,xout_loss(iCO,:),'m--')
plot(O2C,xout_loss(iCO2,:),'c--')
plot(O2C,xout_loss(iO2,:),'y--')
plot([0 1],[0 0],'k')
plot([0 1],[0 0],'k:')
plot([0 1],[0 0],'k--')
legend('H_2','H_2O','CH_4','CO','CO_2','O_2')
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Equilibrium Mole Fraction')
text(.02,0.62,sprintf('Autothermal Reforming'))
text(.02,0.56,'Conditions: 10 bar, 450\circC preheat')
text(.02,0.50,sprintf('Water/Carbon Mole Ratio:  %.1f',H2O_C))
axis([0 1 0 0.5])
hold off
plotfixer

figure(3)
clf
hold on
plot(O2C,CGE,'b')
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('CGE')
axis([0 1 0 120])
hold off
plotfixer

figure(4)
clf
hold on
plot(O2C,Syngas_yield,'b')
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Syngas_yield')
axis([0 1 0 3])
hold off
plotfixer


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






