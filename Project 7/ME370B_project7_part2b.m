
%Create the following graphs: T-s, T-h, P-h, h-s spline fits
% Include condentherm condenbar, inflection point, and critical point

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

T_steps = length(Tspline);

% Generate Secondary Properties (h, s)
h_dome = zeros(1,length(Tspline));
s_dome = zeros(1,length(Tspline));
for i = 1:T_steps
    T = Tspline(i)
    h_dome(i) = h_crT(c, rho_vap_dome_spl(i), T);
    s_dome(i) = s_crT(c, rho_vap_dome_spl(i), T);
end

hcritair = h_crT(c, rcritair, Tcritair);
hinfl = h_crT(c, rinfl, Tinfl);
hmaxcondentherm = h_crT(c, rmaxcondentherm, Tmaxcondentherm);
hmaxcondenbar = h_crT(c, rmaxcondenbar, Tmaxcondenbar);
scritair = s_crT(c, rcritair, Tcritair);
sinfl = s_crT(c, rinfl, Tinfl);
smaxcondentherm = s_crT(c, rmaxcondentherm, Tmaxcondentherm);
smaxcondenbar = s_crT(c, rmaxcondenbar, Tmaxcondenbar);

% Graph 1: T-s
figure(1)
plot(s_dome, Tspline)
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
plot(h_dome, Tspline)
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
semilogy(h_dome, Pspline)
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
plot(s_dome, h_dome)
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
