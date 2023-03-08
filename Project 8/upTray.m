function [vapout, liqin] = upTray(vapin, liqout)
% gives the vapor out and the liquid in given the strutures vapor in and
% liquid out, which have all information about the two streams

% Molar masses
MW_O2 = .031999; % (kg/mol)
MW_N2 = .0280134; % (kg/mol)
MW_Ar = .039948; % (kg/mol)

% Guess mass flow rate of vapout
vapout.mdot = (vapin.mdot + liqout.mdot)/2;

mdotinc = vapout.mdot*1e-3; % increment on mass flow rate
Qtoler = 1;
mdottoler = 1e-3;
mdotlast = vapout.mdot;

liqin.P = vapin.P;
vapout.P = vapin.P;

imax = 15; % max NR iterations
for i = 1:1:imax
    vapout.mdot

    % vapout is on tie line with liqout
    [vapout.T,~,vapout.r,vapout.c] = Fast_Bubble_cP(liqout.c,liqout.P);
    vapout.h = h_crT(vapout.c,vapout.r,vapout.T);
    
    % get liqin properties from mass/ mole balance
    liqin.mdot = liqout.mdot + vapout.mdot - vapin.mdot;
    
    vapout.ndot = vapout.mdot/(MW_N2*vapout.c(N2)+MW_O2*vapout.c(O2)+MW_Ar*vapout.c(Ar));
    
    c_hold(N2) = vapout.c(N2)*vapout.ndot + liqout.c(N2)*liqout.ndot - vapin.c(N2)*vapin.ndot;
    c_hold(O2) = vapout.c(O2)*vapout.ndot + liqout.c(O2)*liqout.ndot - vapin.c(O2)*vapin.ndot;
    c_hold(Ar) = vapout.c(Ar)*vapout.ndot + liqout.c(Ar)*liqout.ndot - vapin.c(Ar)*vapin.ndot;
    liqin.ndot = sum(c_hold);
    liqin.c = c_hold/liqin.ndot;
    
    % check if this is right method for enthalpy??
    [liqin.T,liqin.r,~,~] = Fast_Bubble_cP(liqin.c,liqin.P);
    liqin.h = h_crT(liqin.c,liqin.r,liqin.T);
    
    Q = liqin.h*liqin.mdot + vapin.h*vapin.mdot - vapout.h*vapout.mdot - liqout.h*liqout.mdot;
    
    if ((abs(Q) < Qtoler)&&(abs(vapout.mdot-mdotlast)/vapout.mdot < mdottoler))
        return
    end
    
    mdotlast = vapout.mdot;
    
    % Use Newton-Raphson
    % central difference
    % low
    mdotlow = vapout.mdot - mdotinc;
    [Tlow,rllow,rvlow,clow] = Fast_Bubble_cP(liqout.c,liqout.P);
    vapout_hlow = h_crT(clow,rvlow,Tlow);
    liqinmdotlow = liqout.mdot + mdotlow - vapin.mdot;
    ndotlow = mdotlow/(MW_N2*clow(N2)+MW_O2*clow(O2)+MW_Ar*clow(Ar));
    
    c_holdlow(N2) = clow(N2)*ndotlow + liqout.c(N2)*liqout.ndot - vapin.c(N2)*vapin.ndot;
    c_holdlow(O2) = clow(O2)*ndotlow + liqout.c(O2)*liqout.ndot - vapin.c(O2)*vapin.ndot;
    c_holdlow(Ar) = clow(Ar)*ndotlow + liqout.c(Ar)*liqout.ndot - vapin.c(Ar)*vapin.ndot;
    liqinndotlow = sum(c_holdlow);
    liqinclow = c_holdlow/liqinndotlow;
    
    [liqinTlow,liqinrlow,~,~] = Fast_Bubble_cP(liqinclow,liqin.P);
    liqinlowh = h_crT(liqinclow,liqinrlow,liqinTlow);
    
    Qlow = liqinlowh*liqinmdotlow + vapin.h*vapin.mdot - vapout_hlow*mdotlow - liqout.h*liqout.mdot;
    
    % high
    mdothigh = vapout.mdot + mdotinc;
    [Thigh,rlhigh,rvhigh,chigh] = Fast_Bubble_cP(liqout.c,liqout.P);
    vapout_hhigh = h_crT(chigh,rvhigh,Thigh);
    liqinmdothigh = liqout.mdot + mdothigh - vapin.mdot;
    ndothigh = mdothigh/(MW_N2*chigh(N2)+MW_O2*chigh(O2)+MW_Ar*chigh(Ar));
    
    c_holdhigh(N2) = chigh(N2)*ndothigh + liqout.c(N2)*liqout.ndot - vapin.c(N2)*vapin.ndot;
    c_holdhigh(O2) = chigh(O2)*ndothigh + liqout.c(O2)*liqout.ndot - vapin.c(O2)*vapin.ndot;
    c_holdhigh(Ar) = chigh(Ar)*ndothigh + liqout.c(Ar)*liqout.ndot - vapin.c(Ar)*vapin.ndot;
    liqinndothigh = sum(c_holdhigh);
    liqinchigh = c_holdhigh/liqinndothigh;
    
    [liqinThigh,liqinrhigh,~,~] = Fast_Bubble_cP(liqinchigh,liqin.P);
    liqinhighh = h_crT(liqinchigh,liqinrhigh,liqinThigh);
    
    Qhigh = liqinhighh*liqinmdothigh + vapin.h*vapin.mdot - vapout_hhigh*mdothigh - liqout.h*liqout.mdot;
    
    dQdmdot = (Qhigh-Qlow)/(2*mdotinc);
    
    vapout.mdot = vapout.mdot - Q/dQdmdot;
    
end

disp('Netwon-Raphson failed to converge for Upward Tray');



end