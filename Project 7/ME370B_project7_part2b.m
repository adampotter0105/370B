
%Create the following graphs: T-s, T-h, P-h, h-s spline fits
% Include condentherm condenbar, inflection point, and critical point

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

Tmax = 120;
Tmin = Tlower;
T_steps = 20;
Tlist = linspace(Tmin, Tmax, T_steps);

% Creeate and save data
% P_dew = NaN(1,T_steps);
% P_bub = NaN(1,T_steps);
% rv_dew = NaN(1,T_steps);
% rl_dew = NaN(1,T_steps);
% rl_bub = NaN(1,T_steps);
% rv_bub = NaN(1,T_steps);
% X_dew = NaN(N,T_steps);
% X_bub = NaN(N,T_steps);
% 
% for i = 1:T_steps
%     T = Tlist(i)
%     [P_dew(i), rv_dew(i), rl_dew(i), X_dew(:,i)] = Ideal_Dew_cT(c,T);
%     [P_bub(i), rl_bub(i), rv_bub(i), X_bub(:,i)] = Ideal_Bubble_cT(c,T);
% end
% 
% save 'dew_bubble_data.mat' P_bub P_dew rv_dew rv_bub rl_bub rl_dew X_dew X_bub
load 'dew_bubble_data.mat'

% Generate Secondary Properties (h, s)
h_dew = zeros(1,length(Tlist));
h_bub = zeros(1,length(Tlist));
s_dew = zeros(1,length(Tlist));
s_bub = zeros(1,length(Tlist));
for i = 1:T_steps
    T = Tlist(i)
    h_dew(i) = h_crT(c, rv_dew(i), T);
    s_dew(i) = s_crT(c, rv_dew(i), T);
    h_bub(i) = h_crT(c, rl_bub(i), T);
    s_bub(i) = s_crT(c, rl_bub(i), T);
end

hcritair = h_crT(c, rcritair, Tcritair);
hinfl = h_crT(c, rinfl, Tinfl);
hmaxcondentherm = h_crT(c, rmaxcondentherm, Tmaxcondentherm);
hmaxcondenbar = h_crT(c, rmaxcondenbar, rmaxcondenbar);
scritair = s_crT(c, rcritair, Tcritair);
sinfl = s_crT(c, rinfl, Tinfl);
smaxcondentherm = s_crT(c, rmaxcondentherm, Tmaxcondentherm);
smaxcondenbar = s_crT(c, rmaxcondenbar, rmaxcondenbar);

% Graph 1: T-s
figure(1)
s_q = linspace(min([s_dew,s_bub]), max([s_dew,s_bub]), 100);
plot(s_q, interp1([s_dew,s_bub], [Tlist, Tlist], s_q, 'spline' ))
hold on
scatter(scritair, Tcritair, "k*")
scatter(sinfl, Tinfl, "*r")
scatter(smaxcondenbar, Tmaxcondenbar, "*b")
scatter(smaxcondentherm, Tmaxcondentherm, "*g")
legend(["Spline Fit", "Crit Point", "Inflection Pnt", "Max Condenbar", "Max Condentherm"])
hold off
xlabel("Entropy (J/kg*K)")
ylabel("Temperature (K)")
improvePlot

% Graph 2: T-h
figure(2)
h_q = linspace(min([h_dew,h_bub]), max([h_dew,h_bub]), 100);
plot(h_q, interp1([h_dew,h_bub], [Tlist, Tlist], h_q, 'spline' ))
hold on
scatter(hcritair, Tcritair, "k*")
scatter(hinfl, Tinfl, "*r")
scatter(hmaxcondenbar, Tmaxcondenbar, "*b")
scatter(hmaxcondentherm, Tmaxcondentherm, "*g")
legend(["Spline Fit", "Crit Point", "Inflection Pnt", "Max Condenbar", "Max Condentherm"])
hold off
xlabel("Enthalpy (J/kg)")
ylabel("Temperature (K)")
improvePlot

% Graph 3: P-h
figure(3)
plot(h_q, interp1([h_dew,h_bub], [P_dew, P_bub], h_q, 'spline' ))
hold on
scatter(hcritair, Pcritair, "k*")
scatter(hinfl, Pinfl, "*r")
scatter(hmaxcondenbar, Pmaxcondenbar, "*b")
scatter(hmaxcondentherm, Pmaxcondentherm, "*g")
legend(["Spline Fit", "Crit Point", "Inflection Pnt", "Max Condenbar", "Max Condentherm"])
hold off
xlabel("Enthalpy (J/kg)")
ylabel("Pressure (Pa)")
improvePlot

% Graph 4: h-s
figure(4)
plot(s_q, interp1([s_dew,s_bub], [h_dew, h_bub], s_q, 'spline' ))
hold on
scatter(scritair, hcritair, "k*")
scatter(sinfl, hinfl, "*r")
scatter(smaxcondenbar, hmaxcondenbar, "*b")
scatter(smaxcondentherm, hmaxcondentherm, "*g")
legend(["Spline Fit", "Crit Point", "Inflection Pnt", "Max Condenbar", "Max Condentherm"])
hold off
xlabel("Entropy (J/kg*K)")
ylabel("Enthalpy (J/kg)")
improvePlot
