function [sat_pres, r_sat_liq, r_sat_vap, failed] = Saturation_iT_NR(T)
%Liquid_Spinodal_iT returns the location/density (kg/m^3) of the liquid
%numerical spinodals (outermost locations at which dP/dr goes to zero) for
%any given temperature (K).
%Dongwon Ka

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

% chemical potential min_spec_vol to liquid_spinodal point
% Set limits and step size.
rmax_1 = rupper_i(ispecies);
rmin_1 = Liquid_Spinodal_iT(T);

failed = 0;

% chemical potential vapor_spinodal point to max_spec_vol
rmax_2   = Vapor_Spinodal_iT(T);
rmin_2   = rgtrip_i(ispecies);

% Initialize densities.
r1_guess = (rmax_1-rmin_1)/2;
r2_guess = (rmax_2-rmin_2)/2;

P1 = P_irT(ispecies, r1_guess, T);
P2 = P_irT(ispecies, r2_guess, T);
mu1 = mu_irT(ispecies, r1_guess, T);
mu2 = mu_irT(ispecies, r2_guess, T);
r = [r1_guess; r2_guess]; % parameter vector

% Fitness function
f = @(P1, P2, mu1, mu2) (P1-P2)^2 + (mu1 - mu2)^2;

% Tunable parameters (IF HAVING PROBLEMS INCREASE EPSILON, DECREASE ALPHA)
epsilon = 1500; % Maximum acceptable error (actually small compared to P, mu values)
alpha = 0.05; % Scaler to make steps less agressive (slower)
drmax1 = (rmax_1 - rmin_1)*0.05;
drmax2 = (rmax_2 - rmin_2)*0.05;

% Derivative estimate of f wrt x
dfdr = @(f, r) (f(1.001*r) - f(0.999*r))/(0.002*r); 

it = 0;
while f(P1, P2, mu1, mu2) > epsilon
    it = it + 1;
    %Evaluate at guesses
    r1_guess = r(1);
    r2_guess = r(2);
    P1 = P_irT(ispecies, r1_guess, T);
    P2 = P_irT(ispecies, r2_guess, T);
    mu1 = mu_irT(ispecies, r1_guess, T);
    mu2 = mu_irT(ispecies, r2_guess, T);

    % Pressure and Chemical Potential Error
    dP = P1 - P2; 
    dmu = mu1 - mu2;
    f0 = f(P1, P2, mu1, mu2);

    % Paartial Derivative Estimates
    dP1dr1 = dfdr(@(r) P_irT(ispecies, r, T), r1_guess);
    dmu1dr1 = dfdr(@(r) mu_irT(ispecies, r, T), r1_guess);
    dP2dr2 = dfdr(@(r) P_irT(ispecies, r, T), r2_guess);
    dmu2dr2 = dfdr(@(r) mu_irT(ispecies, r, T), r2_guess);

    % Jacobian and NR Step
    J = [2*dP*dP1dr1 + 2*dmu*dmu1dr1; ...
        -2*dP*dP2dr2 - 2*dmu*dmu2dr2];
    %%
    dr = alpha*f0*1./J;
    if dr(1) > drmax1
        dr(1) = drmax1;
    end

    if dr(2) > drmax2
        dr(2) = drmax2;
    end
    %%

    r = r - dr;

    % Clamp values arund max and min numbers
    if r(1) > rmax_1
        r(1) = rmax_1;
    elseif r(1) < rmin_1
        r(1) = rmin_1;
    end

    if r(2) > rmax_2
        r(2) = rmax_2;
    elseif r(2) < rmin_2
        r(2) = rmin_2;
    end
    %r
    if it>1e4
        fprintf("Saturation_it_NR failed to converge \n")
        failed = 1;
        break
    end
end


%72.1952
    %0.9361
% output
sat_pres  = P_irT(ispecies, r(1),T);
r_sat_liq = r(1);
r_sat_vap = r(2);