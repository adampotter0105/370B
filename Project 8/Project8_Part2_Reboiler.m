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
N2_liqin = zeros(1,length(qualities));
O2_liqin = zeros(1,length(qualities));
Ar_liqin = zeros(1,length(qualities));

for i = 1:length(qualities)
    [Q(i), vapor(i), liquid(i)] = Reboiler_xqP(x_boil, qualities(i), P_boil);
    [vapout_tray(i), liqin_tray(i)] = upTray(vapor(i),liquid(i));
    N2_vapout(i) = vapout_tray(i).c(N2);
    O2_vapout(i) = vapout_tray(i).c(O2);
    Ar_vapout(i) = vapout_tray(i).c(Ar);
    N2_liqin(i) = liqin_tray(i).c(N2);
    O2_liqin(i) = liqin_tray(i).c(O2);
    Ar_liqin(i) = liqin_tray(i).c(Ar);
end

% PLOTS

figure(1)
plot(qualities, N2_vapout, 'r-o');
hold on;
plot(qualities, O2_vapout,'g-o');
plot(qualities, Ar_vapout,'b-o');
plot(qualities, N2_liqin,'r--o');
plot(qualities, O2_liqin,'g--o');
plot(qualities, Ar_liqin,'b--o');
legend('N2 Vapor Out','O2 Vapor Out','Ar Vapor Out','N2 Liquid In','O2 Liquid In','Ar Liquid In');
hold off;
improvePlot

