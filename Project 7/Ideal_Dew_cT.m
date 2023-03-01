function [P_dew, rg, rf, x] = Ideal_Dew_cT(c,T)
% return the ideal dew point given a specified mole fraction vector and
% temperature
addpath 'Fundamental Relation Files'
addpath 'Fundamental Relation Data'
addpath 'Mixture Models'
addpath 'Setup Files'
addpath 'Property Files'
addpath 'Procedure Files'

format compact

% Set up the basic storage and load the FR files and mixture model.
% Set the number of components in mixture (N = 2 for binary, N = 3 for ternary).
N = length(c);
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

% Get the inflection point for this composition.
[Tinfl, rinfl] = Pr_Inflection_c(c);
Pinfl = P_crT(c,rinfl,Tinfl);

    Pmax = Pcritair;
    Pmin = Psolidair;
    x = zeros(N,1);
    it = 0;
    failed = 0;
    while abs(sum(x)-1) > 0.001 || isnan(sum(x))
        it = it + 1;
        % Set Pressure
        P = (Pmax-Pmin)/2 + Pmin;
        % Get rv of mixture
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
        rho = [rl_N2, rl_O2, rl_Ar];
        % Get the complement
        MW_air = x(O2)*MW_O2+x(N2)*MW_N2+x(Ar)*MW_Ar; % (g/mol)
        mass(O2) = x(O2)*MW_air/1000; % (kg)
        mass(N2) = x(N2)*MW_air/1000; % (kg)
        mass(Ar) = x(Ar)*MW_air/1000; % (kg)
        V(O2) = mass(O2)/rho(2); % (m^3)
        V(N2) = mass(N2)/rho(1); % (m^3)
        V(Ar) = mass(Ar)/rho(3); % (m^3)
        rlc = sum(mass)/sum(V); % (kg/m^3)
        
        if isnan(sum(x))
            Pmax = 0.9*(Pmax-Pmin) + Pmin; % Reduce Pmax to find area without NaNs
        elseif sum(x) > 1
            Pmax = P; % Make new P lower
        else
            Pmin = P; % Make new P higher
        end

        
        if it>100
            fprintf("Failed to Find Solution \n")
            failed = 1;
            break
        end
    end
    % Storage for plotting
    if failed == 1
        statement = ['Failed to find ideal dew point for T = ',num2str(T)];
        disp(statement);
        P_dew = 0;
        rg = 0;
        rf = 0;
        x = c;
        return
    end
    P_dew    = P;
    rg  = rv;
    rf  = rlc;
    XN2_dew  = x(N2);
    XO2_dew  = x(O2);
    XAr_dew  = x(Ar);

end