function [vapout, liqin] = upTray(vapin, liqout)
% gives the vapor out and the liquid in given the strutures vapor in and
% liquid out, which have all information about the two streams

% Inputs: vapin: P, mdot, ndot, h, c (1xN array)
%           liqout: P, mdot, ndot, h, c (1xN array)
% Outputs: vapout: P, mdot, ndot, h, c (1xN array)
%           liqin: P, mdot, ndot, h, c (1xN array)

% Molar masses
MW_O2 = .031999; % (kg/mol)
MW_N2 = .0280134; % (kg/mol)
MW_Ar = .039948; % (kg/mol)

% Find mass flow using quality of output
N = length(vapin.c);
if N == 2 % Accomodate binary and ternary air
    MW = [MW_N2, MW_O2];
else 
    MW = [MW_N2, MW_O2, MW_Ar];
end

% Guess mass flow rate of vapout
vapout.mdot = vapin.mdot;

% Bounds on vapout.mdot
vap_mdot_low = 0;
vap_mdot_high = 10*vapout.mdot; % MAY NEED TO CHANGE

% NR Optimization Parameters
mdotinc = vapout.mdot*1e-3; % increment on mass flow rate
Qtoler = 1; % Evaluation criteria
mdotlast = vapout.mdot; % Initialized as the guess

% Column in assumed to be isobaric
liqin.P = vapin.P;
vapout.P = vapin.P;

imax = 15; % max NR iterations

%%  Functions that don't need to be interated on
% vapout is on tie line with liqout
[vapout.T,~,vapout.r,vapout.c] = Fast_Bubble_cP(liqout.c,liqout.P);
vapout.h = h_crT(vapout.c,vapout.r,vapout.T);

%% Begin NR Loop
for i = 1:1:imax
    vapout.mdot
    % Find capout moles from mass
    vapout.ndot = vapout.mdot/dot(vapout.c,MW);
    
    % Mass in = mass out, find liqin mass from guess vapout
    liqin.mdot = liqout.mdot + vapout.mdot - vapin.mdot;
    
    % Moles in = moles out
    liqin.c = vapout.c*vapout.ndot + liqout.c*liqout.ndot - vapin.c*vapin.ndot;
    liqin.c = liqin.c/sum(liqin.c); % Normalize composition
    liqin.ndot = liqin.mdot/dot(liqin.c, MW);
    
    % check if this is right method for enthalpy??
    [liqin.T,liqin.r,~,~] = Fast_Bubble_cP(liqin.c,liqin.P);
    liqin.h = h_crT(liqin.c,liqin.r,liqin.T);
    
    % Find net energy imbalance, for adiabatic tray should be zero
    Q = liqin.h*liqin.mdot + vapin.h*vapin.mdot - vapout.h*vapout.mdot - liqout.h*liqout.mdot;
    
    % Check if successfully converged
    if abs(Q) < Qtoler
        return
    end
    
    mdotlast = vapout.mdot;
    
    % Use Newton-Raphson
    % central difference
    %% LOW INCREMENT
    mdotlow = vapout.mdot - mdotinc;
    % Find capout moles from mass
    vapout_ndot = mdotlow/dot(vapout.c,MW);
    
    % Mass in = mass out, find liqin mass from guess vapout
    liqin_mdot = liqout.mdot + mdotlow - vapin.mdot;
    
    % Moles in = moles out
    liqin_c = vapout.c*vapout_ndot + liqout.c*liqout.ndot - vapin.c*vapin.ndot;
    liqin_c = liqin_c/sum(liqin_c); % Normalize composition

    % check if this is right method for enthalpy??
    [liqin_T,liqin_r,~,~] = Fast_Bubble_cP(liqin_c,liqin.P);
    liqin_h = h_crT(liqin_c,liqin_r,liqin_T);
    
    % Find net energy imbalance, for adiabatic tray should be zero
    Qlow = liqin_h*liqin_mdot + vapin.h*vapin.mdot - vapout.h*mdotlow - liqout.h*liqout.mdot;
    
    %% HIGH INCREMENT
    mdothigh = vapout.mdot + mdotinc;
    % Find capout moles from mass
    vapout_ndot = mdothigh/dot(vapout.c,MW);
    
    % Mass in = mass out, find liqin mass from guess vapout
    liqin_mdot = liqout.mdot + mdothigh - vapin.mdot;
    
    % Moles in = moles out
    liqin_c = vapout.c*vapout_ndot + liqout.c*liqout.ndot - vapin.c*vapin.ndot;
    liqin_c = liqin_c/sum(liqin_c); % Normalize composition

    % check if this is right method for enthalpy??
    [liqin_T,liqin_r,~,~] = Fast_Bubble_cP(liqin_c,liqin.P);
    liqin_h = h_crT(liqin_c,liqin_r,liqin_T);
    
    % Find net energy imbalance, for adiabatic tray should be zero
    Qhigh = liqin_h*liqin_mdot + vapin.h*vapin.mdot - vapout.h*mdothigh - liqout.h*liqout.mdot;


    % Perform increment on parameter
    dQdmdot = (Qhigh-Qlow)/(2*mdotinc);
    
    vapout.mdot = vapout.mdot - Q/dQdmdot;

    % Bisect if out of bounds
    if vapout.mdot < vap_mdot_low
        vapout.mdot = 0.5*(mdotlast - vap_mdot_low) + vap_mdot_low;
    elseif vapout.mdot > vap_mdot_high
        vapout.mdot = 0.5*(vap_mdot_high - mdotlast) + mdotlast;
    end
    
end

disp('Netwon-Raphson failed to converge for Upward Tray');



end