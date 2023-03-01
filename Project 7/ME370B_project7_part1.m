% ME370B: Modeling and Advanced Concepts
% Project 7 - Part 1
% Andy Huynh

% Make a pressure-density diagram.  This file is useful for exploring how the P_crT
% function behaves in various regions (freezing line, P-r inflection point, etc.).
% This version is set up for the O2-N2-AR system and defaults to air.
% C.F. Edwards, 2/19/12 

addpath 'Fundamental Relation Files'
addpath 'Fundamental Relation Data'
addpath 'Mixture Models'
addpath 'Setup Files'
addpath 'Property Files'
addpath 'Procedure Files'

clear all
format compact
fprintf('\n************************************************************\n')

% Set up the basic storage and load the FR files and mixture model.
% Set the number of components in mixture (N = 2 for binary, N = 3 for ternary).
N = 3;
Setup_Air_Props;

% Set fixed-point values that are specific to air.  See Lemmon et al. 2000.
% The composition is: N2:O2:Ar = 0.7812:0.2096:0.0092
% Be sure to change these as needed if you use engineering air (79:21).
% And ignore these if you use an arbitrary composition.  They are for air.
Mair = 28.9586;                 % kg/kmol
Tmaxcondentherm = 132.6312;     % K
Pmaxcondentherm = 3.78502e6;    % Pa
rmaxcondentherm = 10.4477*Mair; % kg/m3
Tmaxcondenbar   = 132.6035;     % K
Pmaxcondenbar   = 3.7891e6;     % Pa
rmaxcondenbar   = 11.0948*Mair; % kg/m3
Tcritair        = 132.5306;     % K
Pcritair        = 3.7860e6;     % Pa
rcritair        = 11.8308*Mair; % kg/m3
% Bottom of dome conditions for air:
Tsolidair       = 59.75;        % K
Psolidair       = 5265;         % Pa
rsolidair       = 33.067*Mair;  % kg/m3
% Upper limit conditions for the model:
Tupper          = 870;          % K
Pupper          = 100e6;        % Pa
rupper          = rsolidair;
% Lower limit to stay just above solid air:
Tlower          = 60;           % K
% Molar masses
MW_O2 = 31.999; % (g/mol)
MW_N2 = 28.0134; % (g/mol)
MW_Ar = 39.948; % (g/mol)

% Set an air composition.
c = zeros(N,1);
if(N == 3)
    % Use ternary air.
    c(O2) = 0.2096;
    c(Ar) = 0.0092;
    c(N2) = 1 - c(O2) - c(Ar);
end
if(N == 2)
    % Use binary air.
    c(O2) = 0.21;
    c(N2) = 1 - c(O2);
end

% Get the inflection point for this composition.
[Tinfl rinfl] = Pr_Inflection_c(c)
Pinfl = P_crT(c,rinfl,Tinfl)

% Set boundaries for plot.
Tmax = Tinfl;
Tmin = Tlower;
steps = 10;
dT = (Tmax - Tmin)/steps;
i = 1;
% Initialize some variables to store data.
Tplot   = zeros(steps+1,1); % (K)
Pplot   = zeros(steps+1,1); % (Pa)
rfiplot = zeros(steps+1,1); % ideal solution bubble point (kg/m^3)
rgiplot = zeros(steps+1,1); % ideal solution dew point (kg/m^3)
rfcplot = zeros(steps+1,1); % complement solution bubble point (kg/m^3)
rgcplot = zeros(steps+1,1); % complement solution dew point (kg/m^3)
XN2plot = zeros(steps+1,1); % mole fraction N2
XO2plot = zeros(steps+1,1); % mole fraction O2
XArplot = zeros(steps+1,1); % mole fraction Ar

for T=Tmin:dT:Tmax
    T
    P = Psolidair;
    x = zeros(N,1);
    while abs(sum(x)-1) > 0.001
        % Update P for next iteration
        P = P + (Pcritair - Psolidair)/10000;
        % Get rv & rl of mixture
        rv = rv_cTP(c,T,P);
        % Get the chemical potentials of each species in mixture c at rv,T
        mu_mix_N2 = mui_icrT(N2,c,rv,T); % (J/kmol-species-i)
        mu_mix_O2 = mui_icrT(O2,c,rv,T);
        mu_mix_Ar = mui_icrT(Ar,c,rv,T);
        % Get the bubble point of each species
        rl_O2 = rl_iTP(O2,T,P);
        rl_N2 = rl_iTP(N2,T,P);
        rl_Ar = rl_iTP(Ar,T,P);
        % Get the chemical potential of each species as pure fluid
        mu_pure_O2 = mu_irT(O2,rl_O2,T); % (J/kmol)
        mu_pure_N2 = mu_irT(N2,rl_N2,T);
        mu_pure_Ar = mu_irT(Ar,rl_Ar,T);
        % Get the mole fractions
        x(O2) = exp((mu_mix_O2 - mu_pure_O2)/(Ru*T));
        x(N2) = exp((mu_mix_N2 - mu_pure_N2)/(Ru*T));
        x(Ar) = exp((mu_mix_Ar - mu_pure_Ar)/(Ru*T));
        % Get the complement
        MW_air = x(O2)*MW_O2+x(N2)*MW_N2+x(Ar)*MW_Ar; % (g/mol)
        mass(O2) = x(O2)*MW_air/1000; % (kg)
        mass(N2) = x(N2)*MW_air/1000; % (kg)
        mass(Ar) = x(Ar)*MW_air/1000; % (kg)
        V(O2) = mass(O2)/rl_O2; % (m^3)
        V(N2) = mass(N2)/rl_N2; % (m^3)
        V(Ar) = mass(Ar)/rl_Ar; % (m^3)
        rlc = sum(mass)/sum(V); % (kg/m^3)
    end
    % Storage for plotting
    Tplot(i)    = T;
    Pplot(i)    = P;
    rgiplot(i)  = rv;
    rlcplot(i)  = rlc;
    XN2plot(i)  = x(N2);
    XO2plot(i)  = x(O2);
    XArplot(i)  = x(Ar);
    i = i+1;
end




% figure(1)
% clf
% hold on
% plot(rinfl,Pinfl/1e6,'kd')
% plot([0 rupper],[Psolidair Psolidair]/1e6,'k--')
% plot(rgsplot,Pgsplot/1e6,'k.',rfsplot,Pfsplot/1e6,'k.')
% plot(rinfl,Pinfl/1e6,'ko')
% %plot([rupper rupper],[-30 110],'r--')
% %plot([0 1000],[100 100],'r--')
% legend('P-\rho Inflection','Solid Air','Location','SouthWest')
% hold off
% xlabel('Density (kg/m^3)')
% ylabel('Pressure (MPa)')
% axis([0 rupper -0.5*Pinfl/1e6 2*Pinfl/1e6])
% if(N == 3)
%     title(sprintf('%.3f N2, %.3f O2, %.3f Ar',c(N2),c(O2),c(Ar)))
% end
% if(N == 2)
%     title(sprintf('%.3f N2, %.3f O2',c(N2),c(O2)))
% end