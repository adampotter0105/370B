% ME370B: Modeling and Advanced Concepts
% Project 7 - Part 2a
% Andy Huynh

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
c(2) = 0.2096; c(3) = 0.0092; c(1) = 1 - c(2) - c(3);
N = length(c);
Setup_Air_Props;

% Set fixed-point values that are specific to air.  See Lemmon et al. 2000.
% The composition is: N2:O2:Ar = 0.7812:0.2096:0.0092
% Be sure to change these as needed if you use engineering air (79:21).
% And ignore these if you use an arbitrary composition.  They are for air.
Mair            = 28.9586;      % kg/kmol
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
% Molar masses:
MW_N2 = 28.0134;    % (g/mol)
MW_O2 = 31.999;     % (g/mol)
MW_Ar = 39.948;     % (g/mol)

% Get the inflection point for this composition.
[Tinfl rinfl]   = Pr_Inflection_c(c);
Pinfl           = P_crT(c,rinfl,Tinfl);

% Get the data.
[T_locus, r_locus, Ti, Pfigc, Pgifc, xfigc, xgifc, rfi, rgi, rfc, rgc] = Vapor_Dome_c(c);

% Plot them figures!
figure(1)
hold on
box on
plot(rmaxcondentherm,Tmaxcondentherm,'kd')
plot(rmaxcondenbar,Tmaxcondenbar,'ksquare')
plot(rcritair,Tcritair,'ko')
plot(rinfl,Tinfl,'rd')
x = [rgi; rfi; rcritair];
y = [Ti; Ti; Tcritair];
X = linspace(min(rgi),max(rfi),1000);
Y = spline(x,y,X);
plot(X,Y,'k')
x = [rgc; rfc; rcritair];
y = [Ti; Ti; Tcritair];
X = linspace(min(rgc),max(rfc),1000);
Y = spline(x,y,X);
plot(X,Y,'b')
plot(rfi,Ti,'ok',rgi,Ti,'ok')
plot(rfc,Ti,'ob',rgc,Ti,'ob')
legend('Maxcondentherm','Maxcondenbar','Critical Point','P-\rho Inflection','Vapor Dome','Complement')
hold off
xlabel('Density (kg/m^3)')
ylabel('Temperature (K)')

figure(2)
hold on
box on
plot(rmaxcondentherm,Pmaxcondentherm/1e6,'kd')
plot(rmaxcondenbar,Pmaxcondenbar/1e6,'ksquare')
plot(rcritair,Pcritair/1e6,'ko')
plot(rinfl,Pinfl/1e6,'rd')
x = [rgi; rfi; rcritair];
y = [Pgifc/1e6; Pfigc/1e6; Pcritair/1e6];
X = linspace(min(rgi),max(rfi),1000);
Y = spline(x,y,X);
plot(X,Y,'k')
x = [rgc; rfc; rcritair];
y = [Pfigc/1e6; Pgifc/1e6; Pcritair/1e6];
X = linspace(min(rgc),max(rfc),1000);
Y = spline(x,y,X);
plot(X,Y,'b')
plot(rfi,Pfigc/1e6,'ok',rgi,Pgifc/1e6,'ok')
plot(rfc,Pgifc/1e6,'ob',rgc,Pfigc/1e6,'ob')
legend('Maxcondentherm','Maxcondenbar','Critical Point','P-\rho Inflection','Vapor Dome','Complement')
hold off
xlabel('Density (kg/m^3)')
ylabel('Pressure (MPa)')

figure(3)
hold on
box on
x = [rgc; rfc];
y = [xfigc(:,1); xgifc(:,1)];
X = linspace(min(rgc),max(rfc),1000);
Y = spline(x,y,X);
plot(X,Y,'b')
x = [rgc; rfc];
y = [xfigc(:,2); xgifc(:,2)];
X = linspace(min(rgc),max(rfc),1000);
Y = spline(x,y,X);
plot(X,Y,'g')
x = [rgc; rfc];
y = [xfigc(:,3)*10; xgifc(:,3)*10];
X = linspace(min(rgc),max(rfc),1000);
Y = spline(x,y,X);
plot(X,Y,'r')
plot(rfc,xgifc(:,1),'ob',rgc,xfigc(:,1),'ob')
plot(rfc,xgifc(:,2),'og',rgc,xfigc(:,2),'og')
plot(rfc,xgifc(:,3)*10,'or',rgc,xfigc(:,3)*10,'or')
legend('N2','O2','10*Ar')
hold off
xlabel('Complementary Phase Density (kg/m^3)')
ylabel('Complementary Phase Mole Fraction')