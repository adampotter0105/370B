% Ideal Dew/Bubble Point with cT input 

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
N     = 3;
Setup_Air_Props;
global Ru;

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

% % Get the inflection point for this composition.
[Tinfl, rinfl] = Pr_Inflection_c(c);
Pinfl = P_crT(c,rinfl,Tinfl); 

Tgiven = linspace(60,120,10);
P_dew = zeros(1,length(Tgiven));
rho_dew = zeros(1,length(Tgiven));
rho_dew_comp = zeros(1,length(Tgiven));
delP = (Pcritair-Psolidair)/1000;

for j = 1:length(Tgiven)

    P0 = Psolidair;
    while true
        rv = rv_cTP(c,Tgiven(j),P0);
        mu_mix_O2 = mui_icrT(O2,c,rv,Tgiven(j));
        mu_mix_Ar = mui_icrT(Ar,c,rv,Tgiven(j));
        mu_mix_N2 = mui_icrT(N2,c,rv,Tgiven(j));

        rl_O2 = rl_iTP(O2,Tgiven(j),P0);
        rl_Ar = rl_iTP(Ar,Tgiven(j),P0);
        rl_N2 = rl_iTP(N2,Tgiven(j),P0);
        mu_pure_O2 = mu_irT(O2,rl_O2,Tgiven(j));
        mu_pure_Ar = mu_irT(Ar,rl_Ar,Tgiven(j));
        mu_pure_N2 = mu_irT(N2,rl_N2,Tgiven(j));

        x_O2 = exp((mu_mix_O2-mu_pure_O2)/Ru/Tgiven(j));
        x_Ar = exp((mu_mix_Ar-mu_pure_Ar)/Ru/Tgiven(j));
        x_N2 = exp((mu_mix_N2-mu_pure_N2)/Ru/Tgiven(j));
        
        x_sum = x_O2+x_Ar+x_N2
        error = abs(1-x_sum);
        if error < 0.001
            break
        end
        P0 = P0+delP;
    end
    P_dew(j) = P0;
    rho_dew(j) = rv;
    rho_dew_comp(j) = x_O2*rl_O2 + x_Ar*rl_Ar + x_N2*rl_N2;
end

figure(1)
% temperture density
scatter(rho_dew,Tgiven);
hold on;
scatter(rho_dew_comp,Tgiven);
xlim([0 1300]);
