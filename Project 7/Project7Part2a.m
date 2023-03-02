g% Project 7 Part 2a

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
[Tinfl, rinfl] = Pr_Inflection_c(c);
Pinfl = P_crT(c,rinfl,Tinfl);


Temps = linspace(60,0.98*Tcritair,10);

for i = 1:length(Temps)
    T = Temps(i);
    
    [p,rg,rf,x] = Dew_cT(c,T);
    P_dew(i) = p;
    rho_dew(i) = rg;
    rho_dew_comp(i) = rf;
    XN2_dew(i) = x(N2);
    XO2_dew(i) = x(O2);
    XAr_dew(i) = x(Ar);
    
    [p,rf,rg,x] = Bubble_cT(c,T);
    P_bub(i) = p;
    rho_bub(i) = rf;
    rho_bub_comp(i) = rg;
    XN2_bub(i) = x(N2);
    XO2_bub(i) = x(O2);
    XAr_bub(i) = x(Ar);
end

rho_vap_dome = [rho_dew rcritair flip(rho_bub)];
rho_vap_dome_spl = linspace(rho_dew(1),rho_bub(1),50);
%rho_vap_spl = [rho_dew(1) rcritair rho_bub(1)];
rho_comp = [rho_bub_comp rcritair flip(rho_dew_comp)];
rho_comp_spl = linspace(rho_bub_comp(1),rho_dew_comp(1),50);
%rho_com_spl = [rho_bub_comp(1) rcritair rho_dew_comp(1)];
Temps_vap_dome = [Temps Tcritair flip(Temps)];
%Temps_vap_spl = [Temps(1) Tcritair Temps(length(Temps))];
P_vap_dome_comp = [P_bub Pcritair flip(P_dew)];
P_vap_dome = flip(P_vap_dome_comp);
Pspline = interp1(rho_vap_dome,P_vap_dome,rho_vap_dome_spl,'spline');
Pspline_comp = interp1(rho_comp,P_vap_dome_comp,rho_comp_spl,'spline');
Tspline = interp1(rho_vap_dome,Temps_vap_dome,rho_vap_dome_spl,'spline');
Tspline_comp = interp1(rho_comp,Temps_vap_dome,rho_comp_spl,'spline');

XN2 = [XN2_bub c(N2) flip(XN2_dew)];
XO2 = [XO2_bub c(O2) flip(XO2_dew)];
XAr = [XAr_bub c(Ar) flip(XAr_dew)];
XN2_spl = interp1(rho_comp,XN2,rho_comp_spl,'spline');
XO2_spl = interp1(rho_comp,XO2,rho_comp_spl,'spline');
XAr_spl = interp1(rho_comp,XAr,rho_comp_spl,'spline');
figure(1)

plot(rcritair,Tcritair,'rd');
hold on;
plot(rmaxcondentherm,Tmaxcondentherm,'kd')
plot(rmaxcondenbar,Tmaxcondenbar,'ksquare')
plot(rinfl,Tinfl,'ro')
plot(rho_vap_dome_spl,Tspline,'k');
plot(rho_comp_spl,Tspline_comp,'b');
plot(rho_dew,Temps,'ko');
plot(rho_dew_comp,Temps,'bo');
plot(rho_bub,Temps,'ko');
plot(rho_bub_comp,Temps,'bo');
hold off;
xlabel('Density (kg/m^3)');
ylabel('Temperature (K)');
title('Ternary Air');
legend('Critical Point','MaxCondentherm','MaxCondenbar','Inflection Point','Vapor Dome','Complement');


figure(2)

plot(rcritair,Pcritair/1e6,'rd');
hold on;
plot(rmaxcondentherm,Pmaxcondentherm/1e6,'kd')
plot(rmaxcondenbar,Pmaxcondenbar/1e6,'ksquare')
plot(rinfl,Pinfl/1e6,'ro')
plot(rho_vap_dome_spl,Pspline/1e6,'k');
plot(rho_comp_spl,Pspline_comp/1e6,'b');
plot(rho_dew,P_dew/1e6,'ko');
plot(rho_dew_comp,P_dew/1e6,'bo');
plot(rho_bub,P_bub/1e6,'ko');
plot(rho_bub_comp,P_bub/1e6,'bo');
hold off;
xlabel('Density (kg/m^3)');
ylabel('Pressure (MPa)');
title('Ternary Air');
legend('Critical Point','MaxCondentherm','MaxCondenbar','Inflection Point','Vapor Dome','Complement');

figure(3)

plot(rho_comp_spl,XN2_spl,'r');
hold on;
plot(rho_comp_spl,XO2_spl,'g');
plot(rho_comp_spl,10*XAr_spl,'b');
scatter(rho_comp, XN2, 'r')
scatter(rho_comp, XO2, 'g')
scatter(rho_comp, 10*XAr, 'b')
hold off;
xlabel('Complimentary Phase Density (kg/m^3)');
ylabel('Complimentary Phase Mole Fraction');
title('Ternary Air');
legend(["Nitrogen", "Oxygen", "10*Argon"]);