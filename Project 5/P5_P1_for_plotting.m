% Gasification using ideal gases in gasification_small.xml
% Does the problem backwards and uses general HCs or oxygenates.
% C.F. Edwards, 1/25/09

clear all
format compact
fprintf('\n************************************************************\n')

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
iNO  = speciesIndex(gas,'NO');
iNO2  = speciesIndex(gas,'NO2');
iN2O  = speciesIndex(gas,'N2O');
iHCN  = speciesIndex(gas,'HCN');
iAR  = speciesIndex(gas,'AR');

iC   = elementIndex(gas,'C');
iH   = elementIndex(gas,'H');
iO   = elementIndex(gas,'O');
iN   = elementIndex(gas,'N');
iS   = elementIndex(gas,'S');
iAr  = elementIndex(gas,'AR');

H2O_C  = 1;  

for i=1:101
    O2C(i) = (101-i)/100;
    if i<71
        O2C_Biomass(i) = (101-i)/100;
    end
end

result_NG      = Natural_Gas();
result_Coal    = Coal();
result_Biomass = Biomass();

figure(1)
clf
hold on
plot(O2C,result_NG.ToutC,'b')
plot(O2C,result_NG.ToutC_loss,'b--')
plot(O2C,result_Coal.ToutC,'r')
plot(O2C,result_Coal.ToutC_loss,'r--')
plot(O2C_Biomass,result_Biomass.ToutC,'g','Color',[0 0.5 0])
plot(O2C_Biomass,result_Biomass.ToutC_loss,'g--','Color',[0 0.5 0])
plot([0 1],[0 0],'k')
plot([0 1],[0 0],'k--')
yline(1350,'--','LineWidth',1)
yline(800,'--','LineWidth',1)
yline(450,'--','LineWidth',1)
xline((14*O2C(44)+5*O2C(44))/19,'--','LineWidth',1)
xline((17.3*O2C(31)+4*O2C(32))/21.3,'--','LineWidth',1)
legend('Natural Gas','','High-Ash Bit. Coal','','Veg Biomass','','No Losses','With Losses')
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Equilibrium Temperature (\circC)')
text(.05,1400,sprintf('Coal-Ash Melt Temperature'))
text(.05,850,sprintf('Biomass-Ash Melt Temperature'))
text(.05,500,sprintf('Reactor Feed Temperature'))
scale = axis;
axis([0 1 0 2000])
hold off
plotfixer

figure(2)
clf
hold on
plot([0 1],[0 0],'k')
plot([0 1],[0 0],'k--')
plot([0 1],[0 0],'k:')

plot(O2C,result_NG.xout_loss(iH2,:),'b')
plot(O2C,result_NG.xout_loss(iH2O,:),'g','Color',[0 0.5 0])
plot(O2C,result_NG.xout_loss(iCH4,:),'r')
plot(O2C,result_NG.xout_loss(iCO,:),'m')
plot(O2C,result_NG.xout_loss(iCO2,:),'c')
plot(O2C,result_NG.xout_loss(iO2,:),'y')

plot(O2C,result_Coal.xout_loss(iH2,:),'b--')
plot(O2C,result_Coal.xout_loss(iH2O,:),'g--','Color',[0 0.5 0])
plot(O2C,result_Coal.xout_loss(iCH4,:),'r--')
plot(O2C,result_Coal.xout_loss(iCO,:),'m--')
plot(O2C,result_Coal.xout_loss(iCO2,:),'c--')
plot(O2C,result_Coal.xout_loss(iO2,:),'y--')

plot(O2C_Biomass,result_Biomass.xout_loss(iH2,:),'b:')
plot(O2C_Biomass,result_Biomass.xout_loss(iH2O,:),'g:','Color',[0 0.5 0])
plot(O2C_Biomass,result_Biomass.xout_loss(iCH4,:),'r:')
plot(O2C_Biomass,result_Biomass.xout_loss(iCO,:),'m:')
plot(O2C_Biomass,result_Biomass.xout_loss(iCO2,:),'c:')
plot(O2C_Biomass,result_Biomass.xout_loss(iO2,:),'y:')

xline((14*O2C(44)+5*O2C(44))/19,'--','LineWidth',1)
xline((17.3*O2C(31)+4*O2C(32))/21.3,':','LineWidth',1)

legend('Nat. Gas','Bit Coal','Veg Bio','H_2','H_2O','CH_4','CO','CO_2','O_2')
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Equilibrium Mole Fraction')
axis([0 1 0 0.5])
hold off
plotfixer


figure(3)
clf
hold on
plot(O2C,result_NG.CGE,'b')
plot(O2C,result_Coal.CGE,'b--')
plot(O2C_Biomass,result_Biomass.CGE,'b:')
xline((14*O2C(44)+5*O2C(44))/19,'--','LineWidth',2)
xline((17.3*O2C(31)+4*O2C(32))/21.3,':','LineWidth',2)
yline(100,'--','LineWidth',2)
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Cold Gas Efficiency (%)')
legend('Natural Gas','High-Ash Bit. Coal','Veg Biomass')
axis([0 1 0 120])
hold off
plotfixer

figure(4)
clf
hold on
plot(O2C,result_NG.Syngas_yield,'b')
plot(O2C,result_Coal.Syngas_yield,'b--')
plot(O2C_Biomass,result_Biomass.Syngas_yield,'b:')
xline((14*O2C(44)+5*O2C(44))/19,'--','LineWidth',2)
xline((17.3*O2C(31)+4*O2C(32))/21.3,':','LineWidth',2)
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Molar Syngas Yield (H_2 + CO)/C')
legend('Natural Gas','High-Ash Bit. Coal','Veg Biomass')
axis([0 1 0 3])
hold off
plotfixer



%figure(3)
clf
hold on

semilogy(O2C,result_NG.xout_loss(iCH4,:),'r');
hold on
semilogy(O2C,result_NG.xout_loss(iO2,:),'g');
hold on
semilogy(O2C,result_NG.xout_loss(iN2,:),'b');
hold on
semilogy(O2C,result_NG.xout_loss(iH2O,:),'c');
hold on
semilogy(O2C,result_NG.xout_loss(iCO,:),'m');
hold on
semilogy(O2C,result_NG.xout_loss(iCO2,:),'y');
hold on
semilogy(O2C,result_NG.xout_loss(iH2,:),'k');
hold on
semilogy(O2C,result_NG.xout_loss(iSO2,:),'color',[0 0.5 0]);
hold on
semilogy(O2C,result_NG.xout_loss(iCOS,:),'color',[0.8500 0.3250 0.0980]);
hold on
semilogy(O2C,result_NG.xout_loss(iNH3,:),'color',[0.9290 0.6940 0.1250]);
hold on
semilogy(O2C,result_NG.xout_loss(iH2S,:),'color',[0.4940 0.1840 0.5560]);
hold on
semilogy(O2C,result_NG.xout_loss(iNO,:),'color',[0.4660 0.6740 0.1880]);
hold on
semilogy(O2C,result_NG.xout_loss(iNO2,:),'color',[0.3010 0.7450 0.9330]);
hold on
semilogy(O2C,result_NG.xout_loss(iN2O,:),'color',[0.6350 0.0780 0.1840]);
hold on
semilogy(O2C,result_NG.xout_loss(iHCN,:),'color',[0 0.4470 0.7410]);
hold on
%semilogy(O2C,result_NG.xout_loss(iAR,:),'color',[0.5 0.5 0.5]);

axis([0 1.1 10^(-6) 1])
legend('CH4','O2','N2','H2O','CO','CO2','H2','SO2','COS','NH3','H2S','NO','NO2','N2O','HCN')
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Equilibrium Mole Fraction')
title('Natural Gas')
hold off
plotfixer


semilogy(O2C,result_Coal.xout_loss(iCH4,:),'r');
hold on
semilogy(O2C,result_Coal.xout_loss(iO2,:),'g');
hold on
semilogy(O2C,result_Coal.xout_loss(iN2,:),'b');
hold on
semilogy(O2C,result_Coal.xout_loss(iH2O,:),'c');
hold on
semilogy(O2C,result_Coal.xout_loss(iCO,:),'m');
hold on
semilogy(O2C,result_Coal.xout_loss(iCO2,:),'y');
hold on
semilogy(O2C,result_Coal.xout_loss(iH2,:),'k');
hold on
semilogy(O2C,result_Coal.xout_loss(iSO2,:),'color',[0 0.5 0]);
hold on
semilogy(O2C,result_Coal.xout_loss(iCOS,:),'color',[0.8500 0.3250 0.0980]);
hold on
semilogy(O2C,result_Coal.xout_loss(iNH3,:),'color',[0.9290 0.6940 0.1250]);
hold on
semilogy(O2C,result_Coal.xout_loss(iH2S,:),'color',[0.4940 0.1840 0.5560]);
hold on
semilogy(O2C,result_Coal.xout_loss(iNO,:),'color',[0.4660 0.6740 0.1880]);
hold on
semilogy(O2C,result_Coal.xout_loss(iNO2,:),'color',[0.3010 0.7450 0.9330]);
hold on
semilogy(O2C,result_Coal.xout_loss(iN2O,:),'color',[0.6350 0.0780 0.1840]);
hold on
semilogy(O2C,result_Coal.xout_loss(iHCN,:),'color',[0 0.4470 0.7410]);
hold on
xline((14*O2C(44)+5*O2C(44))/19,'--','LineWidth',2)
%semilogy(O2C,result_NG.xout_loss(iAR,:),'color',[0.5 0.5 0.5]);

axis([0 1.1 10^(-6) 1])
legend('CH4','O2','N2','H2O','CO','CO2','H2','SO2','COS','NH3','H2S','NO','NO2','N2O','HCN')
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Equilibrium Mole Fraction')
title('High-Ash Bituminous Coal')
hold off
plotfixer


semilogy(O2C_Biomass,result_Biomass.xout_loss(iCH4,:),'r');
hold on
semilogy(O2C_Biomass,result_Biomass.xout_loss(iO2,:),'g');
hold on
semilogy(O2C_Biomass,result_Biomass.xout_loss(iN2,:),'b');
hold on
semilogy(O2C_Biomass,result_Biomass.xout_loss(iH2O,:),'c');
hold on
semilogy(O2C_Biomass,result_Biomass.xout_loss(iCO,:),'m');
hold on
semilogy(O2C_Biomass,result_Biomass.xout_loss(iCO2,:),'y');
hold on
semilogy(O2C_Biomass,result_Biomass.xout_loss(iH2,:),'k');
hold on
semilogy(O2C_Biomass,result_Biomass.xout_loss(iSO2,:),'color',[0 0.5 0]);
hold on
semilogy(O2C_Biomass,result_Biomass.xout_loss(iCOS,:),'color',[0.8500 0.3250 0.0980]);
hold on
semilogy(O2C_Biomass,result_Biomass.xout_loss(iNH3,:),'color',[0.9290 0.6940 0.1250]);
hold on
semilogy(O2C_Biomass,result_Biomass.xout_loss(iH2S,:),'color',[0.4940 0.1840 0.5560]);
hold on
semilogy(O2C_Biomass,result_Biomass.xout_loss(iNO,:),'color',[0.4660 0.6740 0.1880]);
hold on
semilogy(O2C_Biomass,result_Biomass.xout_loss(iNO2,:),'color',[0.3010 0.7450 0.9330]);
hold on
semilogy(O2C_Biomass,result_Biomass.xout_loss(iN2O,:),'color',[0.6350 0.0780 0.1840]);
hold on
semilogy(O2C_Biomass,result_Biomass.xout_loss(iHCN,:),'color',[0 0.4470 0.7410]);
hold on
xline((17.3*O2C(31)+4*O2C(32))/21.3,':','LineWidth',2)
%semilogy(O2C,result_NG.xout_loss(iAR,:),'color',[0.5 0.5 0.5]);

axis([0 1 10^(-6) 1])
legend('CH4','O2','N2','H2O','CO','CO2','H2','SO2','COS','NH3','H2S','NO','NO2','N2O','HCN')
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Equilibrium Mole Fraction')
title('Vegetative Biomass')
hold off
plotfixer


function [result] = Natural_Gas()
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
    iH2S = speciesIndex(gas,'NH3');
    
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
    
    % Natural Gas
    NinCH4  = .907;     % Use one mole of Natural Gas with H2S
    NinH2S  = .0005;
    NinCH4  = NinCH4 - NinH2S;
    NinC2H6 = .036;
    NinC3H8 = .019;
    NinN2   = .018;
    NinCO2  = .010;
    NinO2   = .010;
    
    i = 1;
    for O2_C=1:-0.01:0 % higher to lower
        O2C(i) = O2_C;  % Oxygen to carbon ratio
        
        % number of atoms of 1 mole of Natural Gas
        Nin1_NG   = 1;
        num_C_NG  = Nin1_NG * (1*NinCH4 + 2*NinC2H6 + 3*NinC3H8);
        num_H_NG  = Nin1_NG * (4*NinCH4 + 2*NinH2S  + 6*NinC2H6 + 8*NinC3H8);
        num_N_NG  = Nin1_NG * (2*NinN2);
        num_O_NG  = Nin1_NG * (2*NinCO2 + 2*NinO2);
        num_S_NG  = Nin1_NG * (1*NinH2S);
        
        % number of atoms of 1 mol of C in Natural Gas
        num_C     = num_C_NG/num_C_NG;
        num_H     = num_H_NG/num_C_NG;
        num_N     = num_N_NG/num_C_NG;
        num_O     = num_O_NG/num_C_NG;
        num_S     = num_S_NG/num_C_NG;
    
        % Inlet H2O and Air
        Nin1_O2   = O2_C * num_C;
        Nin1_H2O  = H2O_C * num_C;
        Nin1_N2   = 3.76 * Nin1_O2;
        N2_C      = 3.76 * O2_C;
    
        % inlet mass fraction
        %yin1 = zeros(Nsp,1);
        %Min1_NG  = 1 * MW_NG; % mass of 1 mol of NG
        mor_C    = 0;
        M_mf     = num_C*Ma(iC) + num_H*Ma(iH) + num_S*Ma(iS) + (num_O)*Ma(iO) + num_N*Ma(iN);  % [kg/kmol]
        M_ash    = mor_C*Ma(iO);
        M_maf    = M_mf - M_ash;
        mdot_maf = 1*M_maf; % M_maf = 1*mdot_maf [kg/s]
        mdot_ash = 0;
        mdot_slag= 0;
        ash_maf  = mdot_slag/mdot_maf;
    
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
       
        mdot_syngas   = Nin2(iCO)*Mm(iCO) + Nin2(iH2)*Mm(iH2) + Nin2(iH2O)*Mm(iH2O) + Nin2(iSO2)*Mm(iSO2) + Nin2(iN2)*Mm(iN2) + Nin2(iO2)*Mm(iO2);
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
                          
        Qlhv_maf = 46.10e6;                 % LHV is known, J/kg
        Cp_mf   = 2.00e3;                  % Cp is known, J/kg-K
        Q_loss  = 0.02*Qlhv_maf;              % loss, J/kg 2%
        
        Qlhv = Qlhv_maf;
        Cp   = Cp_mf;
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
        
    
        h_slag_delta = 0;
        h_syngas          = mdot_maf/mdot_syngas * (h_maf + (1+ash_maf)*h_mf_delta + O2_C*h_O2_450/M_maf + N2_C*h_N2_450/M_maf + H2O_C*h_H2O_450/M_maf - (ash_maf)*h_slag_delta); %%%%
    
        % Use the syngas composition and backwards-computed enthalpy to find the
        % equilibrium state.
        set(gas,'T',300,'P',100000,'X',Nin2);   % Pre-equilibrate at low T
        equilibrate(gas,'TP');                  % This can soothe a cranky Cantera.
        set(gas,'H',h_syngas,'P',10*100000);  % Set the correct state.
        equilibrate(gas,'HP');                  % Equilibrate for real.
        ToutC(i)  = temperature(gas) - 273.15;
        xout(:,i) = moleFractions(gas);
        
        h_syngas_loss     = mdot_maf/mdot_syngas * (h_maf + (1+ash_maf)*h_mf_delta + O2_C*h_O2_450/M_maf + N2_C*h_N2_450/M_maf + H2O_C*h_H2O_450/M_maf - (ash_maf)*h_slag_delta - Q_loss); %%%%
        set(gas,'T',300,'P',100000,'X',Nin2);   % Pre-equilibrate at low T
        equilibrate(gas,'TP');                  % This can soothe a cranky Cantera.
        set(gas,'H',h_syngas_loss,'P',10*100000);  % Set the correct state.
        equilibrate(gas,'HP');  
        ToutC_loss(i)  = temperature(gas) - 273.15;
        xout_loss(:,i) = moleFractions(gas);
    
        CGE(i)          = mdot_syngas*LHV_mass(xout_loss(:,i))/(mdot_maf*Qlhv)*100;
    
        % mdot_maf + mdot_O2 + mdot_N2 + mdot_H2O = mdot_Syngas
        % ndot_syngas = mdot_syngas/MW_Syngas
        % ndot_maf = 1
        MW_Syngas       = Mm'*xout_loss(:,i);
        Syngas_yield(i) = mdot_syngas/MW_Syngas*(xout_loss(iH2,i)+xout_loss(iCO,i))/1;
    
        i = i+1;
    end

    result.ToutC = ToutC;
    result.xout = xout;
    result.ToutC_loss = ToutC_loss;
    result.xout_loss = xout_loss;
    result.CGE = CGE;
    result.Syngas_yield = Syngas_yield;

end

function [result] = Coal()
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
    
        mor_C    = num_C - num_O_coal;
        %M_maf     = num_C*Ma(iC) + num_H*Ma(iH) + num_S*Ma(iS) + (num_O-mor_C)*Ma(iO) + num_N*Ma(iN);  % [kg/kmol]
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
        Nin2(iH2O)  = Nin1_H2O;
        Nin2(iSO2)  = 1 * num_S;
        Nin2(iN2)   = Nin1_N2 + 1*0.5*num_N;
        Nin2(iO2)   = Nin1_O2 + 1*0.5*num_O - 0.5*Nin2(iCO) - 1*Nin2(iSO2);
        if(Nin2(iO2) < 0)
            Nin2(iH2O) = Nin2(iH2O) - (-1)*Nin2(iO2)*2;
            Nin2(iH2)  = Nin2(iH2)  + (-1)*Nin2(iO2)*2;
            Nin2(iO2)  = 0;
        end
       
        mdot_syngas   = Nin2(iCO)*Mm(iCO) + Nin2(iH2)*Mm(iH2) + Nin2(iH2O)*Mm(iH2O) + Nin2(iSO2)*Mm(iSO2) + Nin2(iN2)*Mm(iN2) + Nin2(iO2)*Mm(iO2);
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
        
            h_syngas         = mdot_maf/mdot_syngas * (h_maf + (1+ash_maf)*h_mf_delta + O2_C*h_O2_450/M_maf + N2_C*h_N2_450/M_maf + (H2O_C)*h_H2O_450/M_maf - (ash_maf)*h_slag_delta); %%%%
        
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

    result.ToutC = ToutC;
    result.xout = xout;
    result.ToutC_loss = ToutC_loss;
    result.xout_loss = xout_loss;
    result.CGE = CGE;
    result.Syngas_yield = Syngas_yield;

end

function [result] = Biomass()

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
    for O2_C=1:-0.01:0.31 % higher to lower
        O2C(i) = O2_C;  % Oxygen to carbon ratio
        
        % number of atoms from Ultimate analysis
        Nin1_coal   = 1;
        num_C_coal  = Nin1_coal * 0.547/Ma(iC);
        num_H_coal  = Nin1_coal * 0.060/Ma(iH);
        num_O_coal  = Nin1_coal * 0.389/Ma(iO);
        num_N_coal  = Nin1_coal * 0.003/Ma(iN);
        num_S_coal  = Nin1_coal * 0.001/Ma(iS);
        
        % number of atoms of 1 mol of C in Natural Gas
        num_C     = num_C_coal/num_C_coal;
        num_H     = num_H_coal/num_C_coal;
        num_O     = num_O_coal/num_C_coal;
        num_N     = num_N_coal/num_C_coal;
        num_S     = num_S_coal/num_C_coal;
    
        % Proximate analysis
        Fix_C     = 0.15;  % < 15 wt% ar   
        Vol_mater = 0.70;  % > 70 wt% maf
        Water     = 0.20;  % 20   wt% ar
        Ash       = 0.015;  % 1.5  wt% ar
        
        % Inlet H2O and Air
        Nin1_O2   = O2_C * num_C;
        Nin1_H2O  = H2O_C * num_C;
        Nin1_N2   = 3.76 * Nin1_O2;
        N2_C      = 3.76 * O2_C;
    
        mor_C    = num_C - num_O_coal;
        %M_maf     = num_C*Ma(iC) + num_H*Ma(iH) + num_S*Ma(iS) + (num_O-mor_C)*Ma(iO) + num_N*Ma(iN);  % [kg/kmol]
        M_maf     = num_C*Ma(iC) + num_H*Ma(iH) + num_S*Ma(iS) + (num_O)*Ma(iO) + num_N*Ma(iN);  % [kg/kmol]
        
        %M_ash    = mor_C*Ma(iO);  
        mdot_maf = 1*M_maf; % M_maf = 1*mdot_maf [kg/s]
        mdot_ash = mdot_maf*Ash/(Vol_mater); 
        mdot_mf  = mdot_maf*(Vol_mater+Ash)/(Vol_mater);
        mdot_slag= mdot_ash;
        ash_maf  = mdot_slag/mdot_maf;
        mdot_mois= mdot_maf*Water/(Vol_mater);
        ndot_mois= mdot_mois/18;
    
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
       
        mdot_syngas   = Nin2(iCO)*Mm(iCO) + Nin2(iH2)*Mm(iH2) + Nin2(iH2O)*Mm(iH2O) + Nin2(iSO2)*Mm(iSO2) + Nin2(iN2)*Mm(iN2) + Nin2(iO2)*Mm(iO2);
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
                          
        Qlhv_maf = 12.0e6;                 % LHV is known, J/kg
        Cp_mf    = 1.3e3;                  % Cp is known, J/kg-K
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
        
        Tm      = 800; % [oC]
        Cp_slag = 1e3;
        Cp_ash  = 1e3;
        hfg_ash = 240e3;
        
        T_predicted_min=100; % guess
        T_predicted_max=2500;
        while (1)
            T_predicted_avg = (T_predicted_max + T_predicted_min)/2;
            if T_predicted_avg<Tm         %%%%%%%%%% T
                h_slag_delta = Cp_ash*(T_predicted_avg-25) + hfg_ash;  %%% ???
            else
                h_slag_delta = Cp_ash*(Tm-25) + hfg_ash + Cp_slag*(T_predicted_avg-Tm);
            end
        
            h_syngas         = mdot_maf/mdot_syngas * (h_maf + (1+ash_maf)*h_mf_delta + O2_C*h_O2_450/M_maf + N2_C*h_N2_450/M_maf + (H2O_C)*h_H2O_450/M_maf - (ash_maf)*h_slag_delta); %%%%
        
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
    
        T_predicted_min=100; % guess
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
    
    result.ToutC = ToutC;
    result.xout = xout;
    result.ToutC_loss = ToutC_loss;
    result.xout_loss = xout_loss;
    result.CGE = CGE;
    result.Syngas_yield = Syngas_yield;

end


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

