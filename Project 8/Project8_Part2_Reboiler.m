clear all

% Set up Air properties
N = 3;
Setup_Air_Props
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
MW_O2 = 0.031999; % (g/mol)
MW_N2 = 0.0280134; % (g/mol)
MW_Ar = 0.039948; % (g/mol) 

Setup_Air_Props

% Call Function to Solve backward Problem

x_boil(O2) = 0.95; x_boil(N2) = 0.025; x_boil(Ar) = 0.025;
P_boil = 1e5;
qualities = 0:0.1:1;

x_in = NaN(3, length(qualities));
T_boil = NaN(1, length(qualities));
T_in = NaN(1, length(qualities));
Q = NaN(1, length(qualities));

N2_vapout = zeros(1,length(qualities));
O2_vapout = zeros(1,length(qualities));
Ar_vapout = zeros(1,length(qualities));
N2_vapin = zeros(1,length(qualities));
O2_vapin = zeros(1,length(qualities));
Ar_vapin = zeros(1,length(qualities));
N2_liqin = zeros(1,length(qualities));
O2_liqin = zeros(1,length(qualities));
Ar_liqin = zeros(1,length(qualities));
N2_liqout = zeros(1,length(qualities));
O2_liqout = zeros(1,length(qualities));
Ar_liqout = zeros(1,length(qualities));
T_liqin = zeros(1,length(qualities));
T_vapin = zeros(1,length(qualities));
T_liqout = zeros(1,length(qualities));
tray_quality = zeros(1,length(qualities));
tray_quality_mole = zeros(1,length(qualities));

for i = 1:length(qualities)
    [Q(i), vapor(i), liquid(i)] = Reboiler_xqP(x_boil, qualities(i), P_boil);
    [vapout_tray(i), liqin_tray(i)] = upTray(vapor(i),liquid(i));
    N2_vapout(i) = vapout_tray(i).c(N2);
    O2_vapout(i) = vapout_tray(i).c(O2);
    Ar_vapout(i) = vapout_tray(i).c(Ar);
    N2_vapin(i) = vapor(i).c(N2);
    O2_vapin(i) = vapor(i).c(O2);
    Ar_vapin(i) = vapor(i).c(Ar);
    N2_liqin(i) = liqin_tray(i).c(N2);
    O2_liqin(i) = liqin_tray(i).c(O2);
    Ar_liqin(i) = liqin_tray(i).c(Ar);
    N2_liqout(i) = liquid(i).c(N2);
    O2_liqout(i) = liquid(i).c(O2);
    Ar_liqout(i) = liquid(i).c(Ar);
    T_liqin(i) = liqin_tray(i).T;
    T_vapin(i) = vapor(i).T;
    T_liqout(i) = liquid(i).T;
    tray_quality(i) = vapout_tray(i).mdot/(vapout_tray(i).mdot + liquid(i).mdot);
    tray_quality_mole(i) = vapout_tray(i).ndot/(vapout_tray(i).ndot + liquid(i).ndot);
end

% PLOTS

figure(1)
plot(qualities, N2_liqin, 'b');
hold on;
plot(qualities, O2_liqin,'r');
plot(qualities, Ar_liqin*10,'g');
plot(qualities, N2_vapout, '--k','LineWidth',1.5);
plot(qualities, N2_liqin, 'k');
plot(qualities, N2_vapout, 'ok');
plot(qualities, N2_vapin, '+k');

plot(qualities, N2_vapout, '--ob','LineWidth',1.5);
plot(qualities, O2_vapout,'--or','LineWidth',1.5);
plot(qualities, Ar_vapout*10,'--og','LineWidth',1.5);
plot(qualities, N2_liqin,'-ob','LineWidth',1.5);
plot(qualities, O2_liqin,'-or','LineWidth',1.5);
plot(qualities, Ar_liqin*10,'-og','LineWidth',1.5);
plot(qualities, N2_vapin, '--+b','LineWidth',1.5);
plot(qualities, O2_vapin, '--+r','LineWidth',1.5);
plot(qualities, Ar_vapin*10, '--+g','LineWidth',1.5);
plot(qualities, N2_liqout, '-+b','LineWidth',1.5);
plot(qualities, O2_liqout, '-+r','LineWidth',1.5);
plot(qualities, Ar_liqout*10, '-+g','LineWidth',1.5);
legend('N2','O2','Ar*10','Vapor','Liquid','Above','Below');
xlabel('Reboiler Outlet Quality');
ylabel('Mole Fraction')
hold off;

figure(2)
plot(qualities, T_liqin, '-o','LineWidth',1.5);
hold on;
plot(qualities, T_vapin, '-o','LineWidth',1.5);
plot(qualities, T_liqout, '-o','LineWidth',1.5);
legend('Liquid Above','Vapor Below','Saturation Outlet');
xlabel('Reboiler Outlet Quality');
ylabel('Temperature (K)');
hold off;

figure(3)
plot(qualities,tray_quality,'--o','LineWidth',1.5);
hold on;
plot(qualities,tray_quality_mole,'-o','LineWidth',1.5);
xlabel('Reboiler Outlet Quality');
ylabel('Tray Outlet Quality');
legend('Mass Quality', 'Mole Quality');
hold on;
