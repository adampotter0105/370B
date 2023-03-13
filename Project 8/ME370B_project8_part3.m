% ME370B: Modeling and Advanced Concepts
% Project 8 - Part 3
% Andy Huynh

clear all

addpath 'Fundamental Relation Files'
addpath 'Fundamental Relation Data'
addpath 'Mixture Models'
addpath 'Setup Files'
addpath 'Property Files'
addpath 'Procedure Files'

% Set up the basic storage and load the FR files and mixture model.
% Set the number of components in mixture (N = 2 for binary, N = 3 for ternary).
N = 2;
Setup_Air_Props;

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

% Set parameters
n_trays = 10;   % number of trays
qC_mass = 0.6;  % quality of condenser
P       = 1e5;  % 1 bar - isobaric process

% Get the target properties at the feed tray from flash input
T1  = 298.15;
P1  = oneatm;
rv1 = rv_cTP(c,T1,P1);
s1  = s_crT(c,rv1,T1);

T2  = T1;
P2  = 200e5;
rv2 = rv_cTP(c,T2,P2);
s2  = s_crT(c,rv2,T2);

T3  = 130;
P3  = P2;
rl3 = rl_cTP(c,T3,P3);
s3  = s_crT(c,rl3,T3);
h3  = h_crT(c,rl3,T3);

P4  = 1e5;
[T4 q4 V4 y4 x4 rg4 rf4] = Flash_zhP(c,h3,P4);
s4l = s_crT(x4,rf4,T4);
s4v = s_crT(y4,rg4,T4);
s4  = s4v*q4 + s4l*(1-q4);
x4  = q4*y4(N2)+(1-q4)*x4(N2);

% Initialize some variables to store data.
Tplot = zeros(1,n_trays+2);
splot = zeros(1,n_trays+2);
xplot = zeros(1,n_trays+2);

% Get the values at the flash (feed tray)
Tplot(end)  = T4;
splot(end)  = s4;
xplot(end)  = x4;

% Moving from top of condenser to bottom
% Guess the composition of vapor outlet
for n = 0.7:0.001:0.999
    % Get the data at the condenser
    y_out_mole(N2)  = n;
    y_out_mole(O2)  = 1-n;
    [Q_out, T_out, rx_out, ry_out, x_out_mole, T_in, y_in, r_in] = Condenser_yqP(y_out_mole, qC_mass, P);
    Tplot(1)        = T_out;
    sl              = s_crT(x_out_mole,rx_out,T_out);
    sv              = s_crT(y_out_mole,ry_out,T_out);
    splot(1)        = sv*qC_mass + sl*(1-qC_mass);
    xplot(1)        = y_out_mole(N2)*qC_mass + x_out_mole(N2)*(1-qC_mass);

    % Rename to input variables to trays
    T_liq_above = T_out;
    x_liq_above = x_out_mole;
    r_liq_above = rx_out;
    T_vap_above = T_in;
    y_vap_above = y_in;
    r_vap_above = r_in;

    % Now run codes to get flow from/to bottom for n_trays
    for i = 2:n_trays+1
        [T_tray, T_liq_below, T_vap_below, x_liq_below, y_vap_below, r_liq_below, r_vap_below, q_tray_mass, q_tray_mole] = Back_Tray(qC_mass, T_liq_above, x_liq_above, r_liq_above, T_vap_above, y_vap_above, r_vap_above, P);
        Tplot(i)        = T_tray;
        sl              = s_crT(x_liq_below,r_liq_below,T_tray);
        sv              = s_crT(y_vap_above,r_vap_above,T_tray);
        splot(i)        = sv*q_tray_mass + sl*(1-q_tray_mass);
        xplot(i)        = y_vap_above(N2)*q_tray_mass + x_liq_below(N2)*(1-q_tray_mass);
        T_liq_above     = T_liq_below;
        x_liq_above     = x_liq_below;
        r_liq_above     = r_liq_below;
        T_vap_above     = T_vap_below;
        y_vap_above     = y_vap_below;
        r_vap_above     = r_vap_below;
    end
    y_bottray = y_vap_above;
    if abs(y_bottray(N2) - y4(N2)) < 0.02
        disp('Success!')
        break
    end
end


% Plots!
% Data for tie lines
y_tie = zeros(1,n_trays+1);
x_tie = zeros(1,n_trays+1);
T_tie = zeros(1,n_trays+1);
sl_tie = zeros(1,n_trays+1);
sv_tie = zeros(1,n_trays+1);
% Add tie lines for flash process
y_ftie = y4(N2);
x_ftie = x4(N2);
T_ftie = T4;
sl_ftie = s_crT(x4,rf4,T4);
sv_ftie = s_crT(y4,rg4,T4);

% Plot 3D T-s-x
figure(1)
hold on
ylabel('Nitrogen Mole Fraction','rotation',0)
xlabel('Specific Entropy (kJ/kg-K)','rotation',0)
zlabel('Temperature (K)')
title('10 trays, 0.60 condenser quality');
view([-10 40])
axis([2 7 0 1 50 300])
grid on
drawnow

% Plot Generated Data
plot3(splot/1e3,xplot,Tplot,'-or','LineWidth',1.5);
drawnow

% Plot tie Lines
plot3([sl_ftie sv_ftie]/1e3, [x_ftie y_ftie], [T_ftie T_ftie], "b--o")
for i = 1:length(T_tie)
    plot3([sl_tie(i) sv_tie(i)]/1e3, [x_tie(i) y_tie(i)], [T_tie(i) T_tie(i)], "b--o")
end
drawnow

load Tsx_Data
% Tell the user the composition planes to be shown.
xN2_planes = clist;

% Add domes at each plane.
for i=1:1:NCS
    plot3(sdome(i,:)/1e3,xN2(i)*ones(1,length(sdome)),Tdome(i,:),'k-','LineWidth',2)
    drawnow
end

% Add dew and bubble points.
for i=1:1:NCS
    plot3(sdew(i,:)/1e3,xN2(i)*ones(1,NPSSC),Tdew(i,:),'b.','LineWidth',2)
    plot3(sbub(i,:)/1e3,xN2(i)*ones(1,NPSSC),Tbub(i,:),'b.','LineWidth',2)
    drawnow
end

% Put splines through the bubble and dew data in the composition direction.
[row col] = size(Tdew);
for j=1
    plot3(sdewSpline(j,:)/1e3,xN2Spline,TdewSpline(j,:),'b-','LineWidth',2)
    plot3(sbubSpline(j,:)/1e3,xN2Spline,TbubSpline(j,:),'b-','LineWidth',2)
    drawnow
end

% Tell the user the pressure surfaces to be shown.
P_surfaces = Plist;

% Add isobars.
for i=1:1:NCS
    % Do the points below the dome.
    for j=1
        plot3(sisoPlow(:,j,i)/1e3,xN2(i)*ones(length(sdome)/2,1),TisoPlow(:,j,i),'b-','LineWidth',2)
        plot3(sisoPhigh(:,j,i)/1e3,xN2(i)*ones(length(sdome)/2,1),TisoPhigh(:,j,i),'b-','LineWidth',2)
        drawnow
    end
end

% Finished with the dynamic plot.  Close it out.
hold off

% Give Side View
view(90,0)
ylim([0 1])
zlim([76 92])





function [Q_out, T_out, rx_out, ry_out, x_out_mole, T_in, y_in, r_in] = Condenser_yqP(y_out_mole, q_mass, P) 
    N = 2;
    % Molar masses:
    MWi(1) = 28.0135;
    MWi(2) = 31.9988;
    
    % Convert outlet composition from mole ratios to mass ratios
    for i = 1:N
        y_out_mass(i) = y_out_mole(i)*MWi(i);
    end
    y_out_mass  = y_out_mass/sum(y_out_mass);   % mass fraction of out gas (sat. vap)
    MW_y_out    = M_c(y_out_mass);              % molar mass of out gas (kg/kmol)
    
    % Get outlet temperature, product densities, and reflux composition
    % along TB-T2-TT tie line
    [T_out ry_out rx_out x_out_mole] = Fast_Dew_cP(y_out_mole, P);
    
    % Convert reflux compositon mole ratios to mass ratios
    for i = 1:N
        x_out_mass(i) = x_out_mole(i)*MWi(i);
    end
    x_out_mass  = x_out_mass/sum(x_out_mass);   % mass fraction of out liq (sat. liq)
    MW_x_out    = M_c(x_out_mass);              % molar mass of out liquid
    
    % Same as converting q_mass to q_mol for species conservation
    % Get mdot out and in
    mdot_y_out      = q_mass; % (kg/s) relative
    mdot_x_out      = 1 - q_mass;
    mdot_in         = mdot_x_out + mdot_y_out;
    
    % Get ndot out
    ndot_y_out      = mdot_y_out/MW_y_out; % (mol/s) or q_mol
    ndot_x_out      = mdot_x_out/MW_x_out;
    
    % Get the composition of sat. vapor feed in
    y_in = ndot_x_out*x_out_mole + ndot_y_out*y_out_mole;
    y_in = y_in/sum(y_in);
    
    % Get density of inlet
    [T_in, r_in, rf_in, x_in] = Fast_Dew_cP(y_in, P);
    
    % Calculates specific heat transfer
    Q_out = h_crT(y_in, r_in, T_in) - q_mass*h_crT(y_out_mole, ry_out, T_out) - (1-q_mass)*h_crT(x_out_mole, rx_out, T_out);
    Q_out = Q_out/1e3;

end

function [T_tray, T_liq_below, T_vap_below, x_liq_below, y_vap_below, r_liq_below, r_vap_below, q_tray_mass, q_tray_mole] = Back_Tray(q_mass, T_liq_above, x_liq_above, r_liq_above, T_vap_above, y_vap_above, r_vap_above, P)
    N = 2;
    % Molar masses:
    MWi(1) = 28.0135;
    MWi(2) = 31.9988;
    
    [T_liq_below rv r_liq_below x_liq_below] = Fast_Dew_cP(y_vap_above, P); % sat. liq. goes out below the tray
    
    % Rename the known:
    T_tray      = T_vap_above;
    r_sat_vapor = r_vap_above;
    y_sat_vapor = y_vap_above;

    % Mass flow rates (known):
    mdot_liq_above = 1-q_mass;
    mdot_sat_vapor = 1;

    % Mole flow rates (known):
    MM_liq_above = x_liq_above*MWi';
    MM_sat_vapor = y_sat_vapor*MWi';
    ndot_liq_above = mdot_liq_above / MM_liq_above;
    ndot_sat_vapor = mdot_sat_vapor / MM_sat_vapor;
    MM_sat_liquid = x_liq_below*MWi';

    % Bisection method for "mdot_sat_liquid":
    mdot_sat_liquid_min = 0.0001;
    mdot_sat_liquid_max = 10; 
    while (1)
        guess = (mdot_sat_liquid_max + mdot_sat_liquid_min)/2;
        mdot_sat_liquid = guess;
        ndot_sat_liquid = mdot_sat_liquid / MM_sat_liquid; 
        mdot_vap_below  = mdot_sat_vapor - mdot_liq_above + mdot_sat_liquid;
    
        y_vap_below = ndot_sat_vapor*y_sat_vapor - ndot_liq_above*x_liq_above + ndot_sat_liquid*x_liq_below;
        y_vap_below = y_vap_below / sum(y_vap_below);
        
        [T_vap_below, r_vap_below, rl, x] = Fast_Dew_cP(y_vap_below, P);

        % Adiabatic Q_dot=0
        zero = mdot_vap_below * h_crT(y_vap_below, r_vap_below, T_vap_below)...
        - mdot_sat_vapor * h_crT(y_sat_vapor, r_sat_vapor, T_tray)...
        + mdot_liq_above * h_crT(x_liq_above, r_liq_above, T_liq_above)...
        - mdot_sat_liquid * h_crT(x_liq_below, r_liq_below, T_tray); 
    
        q_mass
        zero
        guess

        if abs(zero) < 100
            q_tray_mass = mdot_sat_vapor/(mdot_sat_liquid + mdot_sat_vapor);
            q_tray_mole = ndot_sat_vapor/(ndot_sat_liquid + ndot_sat_vapor);
            break
        elseif zero >= 100
            mdot_sat_liquid_max = guess;
        elseif zero <= -100
            mdot_sat_liquid_min = guess;
        end

    end
end