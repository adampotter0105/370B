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
[Tinfl rinfl]   = Pr_Inflection_c(c);
Pinfl           = P_crT(c,rinfl,Tinfl);

% Set boundaries for plot.
Tmax1   = 120; % use 120 instead of Tinfl because of convergence issues
Tmax2   = 100;
Tmin    = Tlower;
steps   = 20;
dT      = (Tmax1 - Tmin)/steps;
steps2  = length(Tmin:dT:Tmax2);
Pmax    = Pcritair;
Pmin    = Psolidair;
Prange  = linspace(Pmin,Pmax,500);
i       = 1;

% Setting the hyperparameters.
m       = 50;      % number of samples
m_elite = 5;       % number of elites samples

% Initialize some variables to store data.
T1plot  = zeros(steps+1,1); % (K)
T2plot  = zeros(steps2,1);  % (K)
P1plot  = zeros(steps+1,1); % (Pa)
P2plot  = zeros(steps2,1);  % (Pa)
rfiplot = zeros(steps2,1);  % ideal solution bubble point (kg/m^3)
rgiplot = zeros(steps+1,1); % ideal solution dew point (kg/m^3)
rfcplot = zeros(steps+1,1); % complement solution bubble point (kg/m^3)
rgcplot = zeros(steps2,1);  % complement solution dew point (kg/m^3)
XfN2plot = zeros(steps+1,1); % mole fraction N2 on bubble side
XfO2plot = zeros(steps+1,1); % mole fraction O2 on bubble side
XfArplot = zeros(steps+1,1); % mole fraction Ar on bubble side
XgN2plot = zeros(steps2,1); % mole fraction N2 on dew side
XgO2plot = zeros(steps2,1); % mole fraction O2 on dew side
XgArplot = zeros(steps2,1); % mole fraction Ar on dew side

fprintf('\nLet''s get started with the dew problem.\n')

% The dew problem.
for T=Tmin:dT:Tmax1
    T
    % Get the starting point.
    for n = 1:length(Prange)
        P = Prange(n);
        % Get rv of mixture
        rv(n)           = rv_cTP(c,T,P);
        if rv(n) == NaN
            continue
        end
        % Get the chemical potentials of each species in mixture c at rv,T
        mu_mix_N2(n)    = mui_icrT(N2,c,rv(n),T); % (J/kmol-species-i)
        mu_mix_O2(n)    = mui_icrT(O2,c,rv(n),T);
        mu_mix_Ar(n)    = mui_icrT(Ar,c,rv(n),T);
        % Get the bubble point of each species
        rl_N2(n)        = rl_iTP(N2,T,P);
        rl_O2(n)        = rl_iTP(O2,T,P);
        rl_Ar(n)        = rl_iTP(Ar,T,P);
        % Get the chemical potential of each species as pure fluid
        mu_pure_N2(n)   = mu_irT(N2,rl_N2(n),T); % (J/kmol)
        mu_pure_O2(n)   = mu_irT(O2,rl_O2(n),T);
        mu_pure_Ar(n)   = mu_irT(Ar,rl_Ar(n),T);
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
    
    fprintf('\nWe''re making progress!\n')
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
            % Get rv of mixture
            rv(s)           = rv_cTP(c,T,P);
            if rv(s) == NaN
                continue
            end
            % Get the chemical potentials of each species in mixture c at rv,T
            mu_mix_N2(s)    = mui_icrT(N2,c,rv(s),T); % (J/kmol-species-i)
            mu_mix_O2(s)    = mui_icrT(O2,c,rv(s),T);
            mu_mix_Ar(s)    = mui_icrT(Ar,c,rv(s),T);
            % Get the bubble point of each species
            rl_N2(s)        = rl_iTP(N2,T,P);
            rl_O2(s)        = rl_iTP(O2,T,P);
            rl_Ar(s)        = rl_iTP(Ar,T,P);
            % Get the chemical potential of each species as pure fluid
            mu_pure_N2(s)   = mu_irT(N2,rl_N2(s),T); % (J/kmol)
            mu_pure_O2(s)   = mu_irT(O2,rl_O2(s),T);
            mu_pure_Ar(s)   = mu_irT(Ar,rl_Ar(s),T);
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
    Vol(N2)     = mass(N2)/rl_N2(id1);      % (m^3)
    Vol(O2)     = mass(O2)/rl_O2(id1);      % (m^3)
    Vol(Ar)     = mass(Ar)/rl_Ar(id1);      % (m^3)
    rlc_best    = sum(mass)/sum(Vol);       % (kg/m^3)
    % Storage for plotting
    T1plot(i)   = T;
    P1plot(i)   = samples(id1);
    rgiplot(i)  = rv(id1);
    rfcplot(i)  = rlc_best;
    XfN2plot(i) = x(N2,id1);
    XfO2plot(i) = x(O2,id1);
    XfArplot(i) = x(Ar,id1);
    i = i+1;
end

fprintf('\nMoving onto the bubble problem...\n')

% The bubble problem.
i = 1;
for T=Tmin:dT:Tmax2
    T
    % Get the starting point.
    for n = 1:length(Prange)
        P = Prange(n);
        % Get rl of mixture
        rl(n)           = rl_cTP(c,T,P);
        if rl(n) == NaN
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
        if (rv_N2(n) == NaN || rv_O2(n) == NaN || rv_Ar(n) == NaN)
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
    
    fprintf('\nWe''re almost there!\n')
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
            if rl(s) == NaN
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
            if (rv_N2(s) == NaN || rv_O2(s) == NaN || rv_Ar(s) == NaN)
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
    T2plot(i)   = T;
    P2plot(i)   = samples(id1);
    rfiplot(i)  = rl(id1);
    rgcplot(i)  = rvc_best;
    XgN2plot(i) = x(N2,id1);
    XgO2plot(i) = x(O2,id1);
    XgArplot(i) = x(Ar,id1);
    i = i+1;
end

fprintf('\nLet''s effing go! Congratulations!\n')

% Plot them figures!
figure(1)
hold on
plot(rmaxcondentherm,Tmaxcondentherm,'kd')
plot(rmaxcondenbar,Tmaxcondenbar,'ksquare')
plot(rcritair,Tcritair,'ko')
plot(rinfl,Tinfl,'rd')
plot(rgiplot,T1plot,'xr')
plot(rfcplot,T1plot,'xb')
plot(rfiplot,T2plot,'xr',rgcplot,T2plot,'xb')
legend('Maxcondentherm','Maxcondenbar','Critical Point','P-\rho Inflection','Ideal Solution','Ideal Comp.')
hold off
xlabel('Density (kg/m^3)')
ylabel('Temperature (K)')

figure(2)
hold on
plot(rmaxcondentherm,Pmaxcondentherm/1e6,'kd')
plot(rmaxcondenbar,Pmaxcondenbar/1e6,'ksquare')
plot(rcritair,Pcritair/1e6,'ko')
plot(rinfl,Pinfl/1e6,'rd')
plot(rgiplot,P1plot/1e6,'xr')
plot(rfcplot,P1plot/1e6,'xb')
plot(rfiplot,P2plot/1e6,'xr',rgcplot,P2plot/1e6,'xb')
legend('Maxcondentherm','Maxcondenbar','Critical Point','P-\rho Inflection','Ideal Solution','Ideal Comp.')
hold off
xlabel('Density (kg/m^3)')
ylabel('Pressure (MPa)')

figure(3)
hold on
plot(rfcplot,XfN2plot,'o')
plot(rfcplot,XfO2plot,'o')
plot(rfcplot,XfArplot*10,'o')
plot(rgcplot,XgN2plot,'o')
plot(rgcplot,XgO2plot,'o')
plot(rgcplot,XgArplot*10,'o')
legend('N2','O2','10*Ar')
hold off
xlabel('Complementary Phase Density (kg/m^3)')
ylabel('Complementary Phase Mole Fraction')