% Make plots that exercise the ideal dew and bubble point functions with air.
% C.F. Edwards, 2/19/12
clear all 
close all

addpath 'Fundamental Relation Files'
addpath 'Fundamental Relation Data'
addpath 'Mixture Models'
addpath 'Setup Files' 
addpath 'Property Files'
addpath 'Procedure Files'

clear all
format compact
fprintf('\n**************************************************************\n')

% Set up the basic storage and load the FR files and mixture model.
% Set the number of components in mixture (N = 2 for binary, N = 3 for ternary).
N = 3;
Setup_Air_Props;

% Set fixed-point values that are specific to air.  See Lemmon et al. 2000.
% The composition is: N2:O2:Ar = 0.7812:0.2096:0.0092
% Be sure to change these as needed if you use engineering air (79:21).
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

% Set the composition.
cair = zeros(N,1);
if(N == 3)
    % Use Lemmon et al. dry, ternary air.
    cair(O2) = 0.2096;
    cair(Ar) = 0.0092;
    cair(N2) = 1 - cair(O2) - cair(Ar);
end
if(N == 2)
    % Use dry, binary (engineering) air.
    cair(O2) = 0.21;
    cair(N2) = 1 - cair(O2);
end
c = cair

%% Backward : Condenser_yqP
% Feed and pressure
P   = 1e5; % pressure 1 bar

% specified outlet composition
y_A_out_mole(N2) = 0.98; 
y_A_out_mole(O2) = 0.01;
y_A_out_mole(Ar) = 0.01;

q_num = 11;
for i=1:q_num
    i
    q_mass(i) = (i-1)/(q_num-1);
    [Q_out(i), T_out(i), r_reflux(i), r_out(i), x_reflux(i,:), T_in(i), y_in(i,:), r_in(i)] = Condenser_yqP(y_A_out_mole, q_mass(i), P);
    [T_sat(i), T_vap_below(i), x_sat_liq(i,:), y_vap_below(i,:), q_tray(i) q_tray_mole(i)] = Back_Tray(q_mass(i), T_out(i), x_reflux(i,:), r_reflux(i), T_in(i), y_in(i,:), r_in(i), P);
end  

y_sat_vapor = y_in;
T_liq_above = T_out;
x_liq_above = x_reflux;


%% plot
figure(1)
plot(q_mass, T_sat,'ko-')
hold on
plot(q_mass, T_liq_above,'b+--')
plot(q_mass, T_vap_below,'bo-')
legend('Sat. Outlet','Liquid Above','Vapor Below')
xlabel('Condenser Outlet Quality(mass)')
ylabel('Temperature (K)')
plotfixer

figure(2)
plot(q_mass, q_tray, 'bo-')
hold on
plot(q_mass, q_tray_mole, 'ro-')
xlabel('Condenser Outlet Quality(mass)')
ylabel('Stage Outlet Quality')
legend('Mass quality','Mole quality')
plotfixer

figure(3)
% for legend
plot(q_mass, inf*ones(length(q_mass)), 'r-')
hold on
plot(q_mass, inf*ones(length(q_mass)), 'b-')
plot(q_mass, inf*ones(length(q_mass)), 'g-')
plot(q_mass, inf*ones(length(q_mass)), 'k-')
plot(q_mass, inf*ones(length(q_mass)), 'k--')
plot(10,10,'k+')
plot(10,10,'ko')

plot(q_mass, y_vap_below(:,1), 'ro-')
plot(q_mass, y_vap_below(:,2), 'bo-')
plot(q_mass, y_vap_below(:,3)*10, 'go-')

plot(q_mass, x_sat_liq(:,1), 'ro--')
plot(q_mass, x_sat_liq(:,2), 'bo--')
plot(q_mass, x_sat_liq(:,3)*10, 'go--')

plot(q_mass, y_sat_vapor(:,1), 'r+-')
plot(q_mass, y_sat_vapor(:,2), 'b+-')
plot(q_mass, y_sat_vapor(:,3)*10, 'g+-')

plot(q_mass, x_liq_above(:,1), 'r+--')
plot(q_mass, x_liq_above(:,2), 'b+--')
plot(q_mass, x_liq_above(:,3)*10, 'g+--')
xlim([0 1])
ylim([0 1])
xlabel('Condenser Outlet Quality(mass)')
ylabel('Mole Fraction')
legend('Nitrogen','','','','','','','','','','', 'Oxygen','','','','','','','','','','', 'Argon*10','','','','','','','','','','', 'Vapor','','','','','','','','','','','Liquid','','','','','','','','','','','Above','Below')
plotfixer

%% Backward Tray
function [T_sat T_vap_below x_sat_liq y_vap_below q_tray q_tray_mole] = Back_Tray(q_mass, T_out, x_reflux, r_reflux, T_in, y_in, r_in, P)
    % molar mass [kg/kmol]
    M_i(1) = 28.0135;
    M_i(2) = 31.9988;
    M_i(3) = 39.9480;
    
    [T_sat_2 rv r_sat_liquid x_sat_liq] = Fast_Dew_cP(y_in, P); % sat. liq. goes out below the tray
    
    % rename the known
    T_liq_above  = T_out;
    r_liq_above  = r_reflux;
    x_liq_above  = x_reflux;
    
    T_sat = T_in;
    r_sat_vapor  = r_in;
    y_sat_vapor  = y_in;

    % mass flow rates (known)
    mdot_liq_above = 1-q_mass;
    mdot_sat_vapor = 1;

    % mole flow rates (known)
    MM_liq_above = x_liq_above*M_i';
    MM_sat_vapor = y_sat_vapor*M_i';
    ndot_liq_above = mdot_liq_above / MM_liq_above;
    ndot_sat_vapor = mdot_sat_vapor / MM_sat_vapor;
    MM_sat_liquid = x_sat_liq*M_i';

    % Bisection method for "mdot_sat_liquid"
    mdot_sat_liquid_min = 0.0001;
    mdot_sat_liquid_max = 10; 
    while (1)
        guess = (mdot_sat_liquid_max + mdot_sat_liquid_min)/2;
        mdot_sat_liquid = guess;
        ndot_sat_liquid = mdot_sat_liquid / MM_sat_liquid; 
        mdot_vap_below  = mdot_sat_vapor - mdot_liq_above + mdot_sat_liquid;
    
        y_vap_below = ndot_sat_vapor*y_sat_vapor - ndot_liq_above*x_liq_above + ndot_sat_liquid*x_sat_liq;
        y_vap_below = y_vap_below / sum(y_vap_below);
        
        [T_vap_below, r_vap_below, rl, x] = Fast_Dew_cP(y_vap_below, P);

        % Adiabatic Q_dot=0
        zero = mdot_vap_below * h_crT(y_vap_below, r_vap_below, T_vap_below)...
        - mdot_sat_vapor * h_crT(y_sat_vapor, r_sat_vapor, T_sat)...
        + mdot_liq_above * h_crT(x_liq_above, r_liq_above, T_liq_above)...
        - mdot_sat_liquid * h_crT(x_sat_liq, r_sat_liquid, T_sat); 
    
        q_mass
        zero
        guess

        if abs(zero) < 100
            q_tray = mdot_sat_vapor/(mdot_sat_liquid + mdot_sat_vapor);
            q_tray_mole = ndot_sat_vapor/(ndot_sat_liquid + ndot_sat_vapor);
            break
        elseif zero >= 100
            mdot_sat_liquid_max = guess;
        elseif zero <= -100
            mdot_sat_liquid_min = guess;
        end

    end
end

%% Condenser_yqP
% Q_out : Heat loss at the top, T_out: temperature of product gas
% r_y_A_out : density of product gas, r_x_A_reflux : density of reflux liq
% x_A_reflux_mole : molar composition of reflux liquid
% T_in : input saturated vapor, y_in : comp. of saturated vapor (input),
% r_in : density of saturated vaopr (input)

function [Q_out, T_out, r_x_A_reflux, r_y_A_out, x_A_reflux_mole, T_in, y_A_in, r_in] = Condenser_yqP(y_A_out_mole, q_mass, P) 
    % molar mass [kg/kmol]
    M_i(1) = 28.0135;
    M_i(2) = 31.9988;
    M_i(3) = 39.9480;
    
    % specified outlet mass ratio
    for i=1:3
        y_A_out_mass(i) = y_A_out_mole(i)*M_i(i); % [kg]
    end
    y_A_out_mass = y_A_out_mass/sum(y_A_out_mass); % mass fraction of out gas (sat. vap)
    y_A_out_MM = M_c(y_A_out_mass);  % molar mass of out gas [kg/kmol]
    
    % outlet temperature, product density, reflux density, and reflux composition 
    %[T_out r_x_A_reflux r_y_A_out x_A_reflux_mole] = Fast_Dew_cP(y_A_out_mole, P);
    [T_out r_y_A_out r_x_A_reflux x_A_reflux_mole] = Fast_Dew_cP(y_A_out_mole, P);

    % reflux compositon
    for i=1:3
        x_A_reflux_mass(i) = x_A_reflux_mole(i)*M_i(i);
    end
    x_A_reflux_mass = x_A_reflux_mass/sum(x_A_reflux_mass); % mass fraction of out liq (sat. liq)
    x_A_reflux_MM = M_c(x_A_reflux_mass);  % molar mass of out liquid
    
    % mass conservation
    mdot_y_A_out    = q_mass; % [kg/s]
    mdot_x_A_reflux = 1 - q_mass;
    mdot_in = mdot_x_A_reflux + mdot_y_A_out;
    
    ndot_y_A_out    = mdot_y_A_out / y_A_out_MM;  % [kmol/s]
    ndot_x_A_reflux = mdot_x_A_reflux / x_A_reflux_MM;
    
    % The composition of sat. vapor feed in
    y_A_in = ndot_x_A_reflux*x_A_reflux_mole + ndot_y_A_out*y_A_out_mole;
    y_A_in = y_A_in/sum(y_A_in);
    
    % density of inlet
    [T_in, r_in, rf_in, x_in] = Fast_Dew_cP(y_A_in, P);
    
    % specific heat calculation
    Q_out = h_crT(y_A_in, r_in, T_in) - q_mass*h_crT(y_A_out_mole, r_y_A_out, T_out) - (1-q_mass)*h_crT(x_A_reflux_mole, r_x_A_reflux, T_out);
    Q_out = Q_out/1e3;

end