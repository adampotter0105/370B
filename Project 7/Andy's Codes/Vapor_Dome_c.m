function [T_locus, r_locus, Ti, rfi, rgi, rfc, rgc] = Vapor_Dome_c(c)
%Vapor_Dome_c returns the T,r locus of the vapor dome, T vector, bubble/dew
%points vectors, and complementary phase bubble/dew vectors
%   c = composition
%   The below is for ternary air:
%   c(2) = 0.2096; c(3) = 0.0092; c(1) = 1 - c(2) - c(3);

% ME370B: Modeling and Advanced Concepts
% Project 7 - Part 2a
% Andy Huynh

addpath 'Fundamental Relation Files'
addpath 'Fundamental Relation Data'
addpath 'Mixture Models'
addpath 'Setup Files'
addpath 'Property Files'
addpath 'Procedure Files'

% Set up the basic storage and load the FR files and mixture model.
% Set the number of components in mixture (N = 2 for binary, N = 3 for ternary).
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
MW_N2           = 28.0134;      % g/mol
MW_O2           = 31.999;       % g/mol
MW_Ar           = 39.948;       % g/mol

% Get the inflection point for this composition.
[Tinfl rinfl]   = Pr_Inflection_c(c);
Pinfl           = P_crT(c,rinfl,Tinfl);

% Setting the boundaries.
Tmax    = 0.98*Tinfl;
Tmin    = Tlower;
steps   = 15;
dT      = (Tmax - Tmin)/steps;

% Intialize some variables to store data.
Ti      = zeros(steps+1,1);
rfi     = zeros(steps+1,1);
rgi     = zeros(steps+1,1);
i       = 1;

% Find the bubble and dew points.
for T=Tmin:dT:Tmax
    T
    [P rf rg y] = Bubble_cT(c,T);
    rfi(i)      = rf;
    rgc(i)      = rg;  
    [P rg rf x] = Dew_cT(c,T);
    rgi(i)      = rg;
    rfc(i)      = rf;
    Ti(i)       = T;
    i = i+1;
end

x = [rgi; rfi];
y = [Ti; Ti];
if N == 3
    x = [x; rcritair];
    y = [y; Tcritair];
end
X = linspace(max(rgi),min(rfi),10000);
Y = spline(x,y,X);
[T_locus idx_max] = max(Y);
r_locus = X(idx_max);

end