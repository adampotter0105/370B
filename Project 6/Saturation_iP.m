function [T_sat,r_liq,r_vap] = Saturation_iP(P)
%Saturation_iP returns the saturation temperature (K), saturated liquid
%density, and saturated vapor density (kg/m^3) for given P (Pa) between the
%triple temperature and near-critical temperature.

% Provide access to support files via the Matlab path.
addpath 'Fundamental_Relation_Files'
addpath 'Fundamental_Relation_Data'
addpath 'Setup_Files'
addpath 'Property_Files'

% Set up the basic storage and load the FR files.
Setup_Props_i;

% Set which of the loaded species you want to work with.  You might want to
% change this around to see what the plots look like for other species.
ispecies = nH2;

Tmax = 21; % approximate max in range of T to evaluate
Tmin = 17; % approximate min in range of T to evaluate
N = 8; % number of data points between Tmax and Tmin

% Iterate a bunch of temperatures to get saturated liquid densities
% (kg/m^3) using saturation function from Part 2b.
i = 1;
for T1 = linspace(Tmin,Tmax,N)
    T(i) = T1;
    [P_sat1(i),r_liq1(i),r_vap1(i)] = Saturation_iT_lookup(T1);
    i = i+1;
end

% Iterate a bunch of temperatures again to get saturated liquid densities
% (kg/m^3) using provided rl_iTP function.
j = 1;
for T2 = linspace(Tmin,Tmax,N)
    [r_liq2(j)] = rl_iTP(ispecies,T2,P);
    j = j+1;
end

% Find minimal point between the two liquid densities.
r_liq_diff = abs((r_liq2-r_liq1));
[diff_min,idx] = min(r_liq_diff);

T_sat = T(idx);
r_liq = rl_iTP(ispecies,T_sat,P);
r_vap = rv_iTP(ispecies,T_sat,P);

end