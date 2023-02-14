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

%     Min1_O2  = Nin1_O2 * Mm(iO2); 
%     Min1_H2O = Nin1_H2O * Mm(iH2O);
%     Min1_N2  = Nin1_N2 * Mm(iN2);
%     Min1_sum = 1*M_mf + Min1_O2 + Min1_H2O + Min1_N2;
%     yin1_mf  = M_mf/Min1_sum;
%     yin1_O2  = Min1_O2/Min1_sum;
%     yin1_H2O = Min1_H2O/Min1_sum;
%     yin1_N2  = Min1_N2/Min1_sum;

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

%     mdot_syngas = Nin(iCO)*Mm(iCO) + Nin(iH2)*Mm(iH2);
%     Min2_syngas = Nin(iCO)*Mm(iCO) + Nin(iH2)*Mm(iH2); 
%     Min2_N2     = Nin(iN2)*Mm(iN2);
%     Min2_O2     = Nin(iO2)*Mm(iO2);
%     Min2_H2O    = Nin(iH2O)*Mm(iH2O);
%     Min2_sum    = Min2_syngas + Min2_H2O + Min2_N2 + Min2_O2;
%     yin2_syngas = Min2_syngas/Min2_sum;
%     yin2_N2     = Min2_N2/Min2_sum;
%     yin2_O2     = Min2_O2/Min2_sum;
%     yin2_H2O    = Min2_H2O/Min2_sum;
    
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
                      
    Qlhv_NG = 46.10e6;                 % LHV is known, J/kg
    Cp_NG   = 2.00e3;                  % Cp is known, J/kg-K
    Q_loss  = 0.02*Qlhv_NG;                       % loss, J/kg 2% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
    Qlhv = Qlhv_NG;
    Cp   = Cp_NG;
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
    h_syngas          = mdot_maf/mdot_syngas * (h_maf + (1+ash_maf)*h_mf_delta + O2_C*h_O2_450/M_maf + N2_C*O2_C*h_N2_450/M_maf + H2O_C*h_H2O_450/M_maf - (ash_maf)*h_slag_delta); %%%%

    % Use the syngas composition and backwards-computed enthalpy to find the
    % equilibrium state.
    set(gas,'T',300,'P',100000,'X',Nin2);   % Pre-equilibrate at low T
    equilibrate(gas,'TP');                  % This can soothe a cranky Cantera.
    set(gas,'H',h_syngas,'P',10*100000);  % Set the correct state.
    equilibrate(gas,'HP');                  % Equilibrate for real.
    
    ToutC(i)  = temperature(gas) - 273.15;
    xout(:,i) = moleFractions(gas);
    
    h_syngas_loss     = mdot_maf/mdot_syngas * (h_maf + (1+ash_maf)*h_mf_delta + O2_C*h_O2_450/M_maf + N2_C*O2_C*h_N2_450/M_maf + H2O_C*h_H2O_450/M_maf - (ash_maf)*h_slag_delta - Q_loss); %%%%
    set(gas,'T',300,'P',100000,'X',Nin2);   % Pre-equilibrate at low T
    equilibrate(gas,'TP');                  % This can soothe a cranky Cantera.
    set(gas,'H',h_syngas_loss,'P',10*100000);  % Set the correct state.
    equilibrate(gas,'HP');  
    ToutC_loss(i)  = temperature(gas) - 273.15;
    xout_loss(:,i) = moleFractions(gas);

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
plot(O2C,xout(iH2,:),'b')
plot(O2C,xout(iH2O,:),'g')
plot(O2C,xout(iCH4,:),'r')
plot(O2C,xout(iCO,:),'m')
plot(O2C,xout(iCO2,:),'c')
plot(O2C,xout(iO2,:),'y')
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