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


for i=1:11
    i
    q_mass(i) = (i-1)/10;
    [Q_out(i), T_out(i), r_reflux(i), r_out(i), x_reflux(i,:), T_in(i), y_in(i,:), r_in(i)] = Condenser_yqP(y_A_out_mole, q_mass(i), P);
end 


%% Plot
figure(1)
plot(q_mass, Q_out,'ro-')
ylim([0 250])
xlabel('Condenser Outlet Quality(mass)')
ylabel('Feed-Specific Heat Rate (kJ/kg)')
plotfixer

figure(2)
plot(q_mass, T_out,'k+--')
hold on
plot(q_mass, T_in,'bo-')
hold off
ylim([77.5 79.0])
xlabel('Condenser Outlet Quality(mass)')
ylabel('Temperature (K)')
legend('Outlet','Inlet')
plotfixer

figure(3)
% inlet feed
plot(q_mass, y_in(:,1),'ro-')
hold on
plot(q_mass, y_in(:,2),'bo-')
plot(q_mass, y_in(:,3)*10,'go-')
% plot for legend
plot(q_mass, inf*ones(11),'ko-')
plot(q_mass, inf*ones(11),'k--')
plot(q_mass, inf*ones(11),'k+--')
% product
plot(q_mass, y_A_out_mole(1)*ones(11),'r--')
plot(q_mass, y_A_out_mole(2)*ones(11),'b--')
plot(q_mass, y_A_out_mole(3)*ones(11)*10,'g--')
% reflux
plot(q_mass, x_reflux(:,1),'r+--')
plot(q_mass, x_reflux(:,2),'b+--')
plot(q_mass, x_reflux(:,3)*10,'g+--')
legend('Nitrogen','Oxygen','Argon*10','Inlet','','','','','','','','','','','','Product','','','','','','','','','','','Reflux')
xlabel('Condenser Outlet Quality(mass)')
ylabel('Mole Fraction')
plotfixer


%% Condenser_yqP
function [Q_out, T_out, r_x_A_reflux, r_y_A_out, x_A_reflux_mole, T_in, y_A_in, r_in] = Condenser_yqP(y_A_out_mole, q_mass, P) 
% molar mass
M_i(1) = 28.0135;
M_i(2) = 31.9988;
M_i(3) = 39.9480;

% specified outlet mass ratio
for i=1:3
    y_A_out_mass(i) = y_A_out_mole(i)*M_i(i);
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
mdot_y_A_out    = q_mass; 
mdot_x_A_reflux = 1 - q_mass;
mdot_in = mdot_x_A_reflux + mdot_y_A_out;

ndot_y_A_out    = mdot_y_A_out / y_A_out_MM;
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





% % guess (input feed composition)
% steps = 200;
% k=1;
% for i=1:steps
%     for j=1:steps-i-1
%         c(k,N2)  = i/steps;
%         c(k,O2)  = j/steps;
%         c(k,Ar)  = 1-c(k,N2)-c(k,O2);
%         k = k+1;
%     end
% end
% 
% for k=1:length(c(:,1))
%     if c(k,1) < x_A_reflux_mole(1)
%         c(k,:)=[0, 0, 0];
%     elseif c(k,1) > y_A_out_mole(1)
%         c(k,:)=[0, 0, 0];
%     end
%     if c(k,2) > x_A_reflux_mole(2)
%         c(k,:)=[0, 0, 0];
%     elseif c(k,2) < y_A_out_mole(2)
%         c(k,:)=[0, 0, 0];
%     end
%     if c(k,3) > x_A_reflux_mole(3)
%         c(k,:)=[0, 0, 0];
%     elseif c(k,3) < y_A_out_mole(3)
%         c(k,:)=[0 0 0];
%     end
% end
% 
% i=1;
% for k=1:length(c(:,1))
%     if c(k,1) ~= 0
%         z(i,:) = c(k,:); 
%         i = i+1;
%     end
% end
% 
% test = 100;
% for i=1:length(z(:,1))
%     i
%     [T_out_2(i) x_A y_A_out_mole_2 rf rg] = Flash_zqP(z(i,:),q_mass,P);
%     if abs(T_out_2(i)-T_out) < test
%         answer = i; 
%         test = abs(T_out_2(i)-T_out);
%     end
% end


% for i=1:11
%     mole_liq(i) = (1-q_mass)/%x_A_out_MM  (sat. liq at 1 bar)
%     mole_vap(i) = q_mass/y_A_out_MM; % sat vapor
% end
% mass_frac_vap_out = q_mass;
% mass_frac_liq_out = 1-q_mass;
% mole_frac_vap_out = q_mass/y_A_out_MM;      % sat vapor
% mole_frac_liq_out = (1-q_mass)/x_A_reflux_MM;  %(sat. liq at 1 bar)
% 
% L_A_reflux = L_A_out*(mole_frac_liq_out/mole_frac_vap_out);




   







