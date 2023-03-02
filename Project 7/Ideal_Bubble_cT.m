function [P_bub, rf, rg, y] = Ideal_Bubble_cT(c,T)
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

Tmax1   = 120; % use 120 instead of Tinfl because of convergence issues
Tmax2   = 100;
Tmin    = Tlower;
steps   = 30;
dT      = (Tmax1 - Tmin)/steps;
steps2  = length(Tmin:dT:Tmax2);
Pmax    = Pcritair;
Pmin    = Psolidair;
Prange  = linspace(Pmin,Pmax,500);

if T > Tmax2
    P_bub = 0;
    rf = 0;
    rg = 0;
    y = c;
    return
end

% Setting the hyperparameters.
m       = 50;      % number of samples
m_elite = 5;       % number of elites samples

% Get the starting point.
for n = 1:length(Prange)
    P = Prange(n);
    % Get rl of mixture
    rl(n)           = rl_cTP(c,T,P);
    if isnan(rl(n))
        continue
    end
    % Get the chemical potentials of each species in mixture c at rl,T
    mu_mix_N2(n)    = mui_icrT(N2,c,rl(n),T); % (J/kmol-species-i)
    mu_mix_O2(n)    = mui_icrT(O2,c,rl(n),T);
    mu_mix_Ar(n)    = mui_icrT(Ar,c,rl(n),T);
    % Get the dew point of each species
    rv_N2(n)        = rv_iTP(N2,T,P);
    rv_O2(n)        = rv_iTP(O2,T,P);
    rv_Ar(n)        = rv_iTP(Ar,T,P);
    if (isnan(rv_N2(n)) || isnan(rv_O2(n)) || isnan(rv_Ar(n)))
        continue
    end
    % Get the chemical potential of each species as pure fluid
    mu_pure_N2(n)   = mu_irT(N2,rv_N2(n),T); % (J/kmol)
    mu_pure_O2(n)   = mu_irT(O2,rv_O2(n),T);
    mu_pure_Ar(n)   = mu_irT(Ar,rv_Ar(n),T);
    % Get the mole fractions
    x(N2,n)         = exp((mu_mix_N2(n) - mu_pure_N2(n))/(Ru*T));
    x(O2,n)         = exp((mu_mix_O2(n) - mu_pure_O2(n))/(Ru*T));
    x(Ar,n)         = exp((mu_mix_Ar(n) - mu_pure_Ar(n))/(Ru*T));
    % Get the sums of x
    sumx(n)         = sum(x(:,n));            
end
[out,idx]   = sort(abs(sumx-1));    % get index of sorted samples
elite_index = idx(1:m_elite);       % retains the best m_elite samples' indices
elites      = Prange(elite_index);  % lists elite samples
avg         = mean(elites);         % calculates new mean of elite samples
covar       = cov(elites);          % calculates new covariance matrix of elite samples

fprintf('\nWe are almost there!\n')
sumx_best = 0;
sumx = [];
while abs(sumx_best-1) > 0.001 % set tolerance
    samples = mvnrnd(avg,covar,m); % create multivariate normal distribution
%         for k = m:-1:1
%             if (samples(k) > Pcritair || samples(k) < Psolidair)
%                 samples(k) = []; % discard infeasible samples
%             end
%         end
    for s = 1:length(samples)
        P = samples(s);
        % Get rl of mixture
        rl(s)           = rl_cTP(c,T,P);
        if isnan(rl(s))
            continue
        end
        % Get the chemical potentials of each species in mixture c at rl,T
        mu_mix_N2(s)    = mui_icrT(N2,c,rl(s),T); % (J/kmol-species-i)
        mu_mix_O2(s)    = mui_icrT(O2,c,rl(s),T);
        mu_mix_Ar(s)    = mui_icrT(Ar,c,rl(s),T);
        % Get the dew point of each species
        rv_N2(s)        = rv_iTP(N2,T,P);
        rv_O2(s)        = rv_iTP(O2,T,P);
        rv_Ar(s)        = rv_iTP(Ar,T,P);
        if (isnan(rv_N2(s)) || isnan(rv_O2(s)) || isnan(rv_Ar(s)))
            continue
        end
        % Get the chemical potential of each species as pure fluid
        mu_pure_N2(s)   = mu_irT(N2,rv_N2(s),T); % (J/kmol)
        mu_pure_O2(s)   = mu_irT(O2,rv_O2(s),T);
        mu_pure_Ar(s)   = mu_irT(Ar,rv_Ar(s),T);
        % Get the mole fractions
        x(N2,s)         = exp((mu_mix_N2(s) - mu_pure_N2(s))/(Ru*T));
        x(O2,s)         = exp((mu_mix_O2(s) - mu_pure_O2(s))/(Ru*T));
        x(Ar,s)         = exp((mu_mix_Ar(s) - mu_pure_Ar(s))/(Ru*T));
        % Get the sums of x
        sumx(s)         = sum(x(:,s));            
    end
    [out,idx]   = sort(abs(sumx-1));    % get index of sorted samples
    elite_index = idx(1:m_elite);       % retains the best m_elite samples' indices
    elites      = samples(elite_index); % lists elite samples
    avg         = mean(elites);         % calculates new mean of elite samples
    covar       = cov(elites);          % calculates new covariance matrix of elite samples
    id1         = idx(1);               % index of best sample
    sumx_best   = sumx(id1);            % closest total mole fractions to 1 so far
end
% Get the complement
MW_air      = x(O2,id1)*MW_O2+x(N2,id1)*MW_N2+x(Ar,id1)*MW_Ar; % (g/mol)
mass(N2)    = x(N2,id1)*MW_air/1000;    % (kg)
mass(O2)    = x(O2,id1)*MW_air/1000;    % (kg)
mass(Ar)    = x(Ar,id1)*MW_air/1000;    % (kg)
Vol(N2)     = mass(N2)/rv_N2(id1);      % (m^3)
Vol(O2)     = mass(O2)/rv_O2(id1);      % (m^3)
Vol(Ar)     = mass(Ar)/rv_Ar(id1);      % (m^3)
rvc_best    = sum(mass)/sum(Vol);       % (kg/m^3)
% Storage for plotting

P_bub   = samples(id1);
rf  = rl(id1);
rg  = rvc_best;
XgN2 = x(N2,id1);
XgO2 = x(O2,id1);
XgAr = x(Ar,id1);
y = [XgN2
     XgO2
     XgAr];

end