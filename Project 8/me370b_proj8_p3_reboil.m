clear all

%% Set up Air properties
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
MW_O2 = 31.999; % (g/mol)
MW_N2 = 28.0134; % (g/mol)
MW_Ar = 39.948; % (g/mol)

Setup_Air_Props

% Functions for saturated liquids and vapors
% [T rl rv y] = Fast_Bubble_cP(x,P,varargin)
% [T rv rl x] = Fast_Dew_cP(y,P,varargin)
%         Tstart  = varargin{1};
%         xstart  = varargin{2};
%         rlstart = varargin{3};
%         rvstart = varargin{4};
% Functions for multi-component flash processes under the vapor dome
% [T x y rf rg] = Flash_zqP(z,q,P,varargin)
% [T q V y x rg rf] = Flash_zhP(z,h,P,varargin)
% [q V y x rg rf] = Flash_zTP(z,T,P,varargin)
%         T  = varargin{1};
%         y  = varargin{2};
%         rf = varargin{3};
%         rg = varargin{4};

% Sate 1: Air at ambient conditions
c(N2) = 0.79; c(O2) = 0.21; % Binary Air
T1 = 298.15;
P1 = oneatm;
rv1 = rv_cTP(c,T1, P1);
s1 = s_crT(c, rv1, T1);

% State 2: Air compressed isothermally to 200 bar
T2 = T1;
P2 = 200e5;
rv2 = rv_cTP(c,T2, P2);
s2 = s_crT(c, rv2, T2);

% State 3: Air Cooled isobarically to 130k
T3 = 130;
P3 = P2;
rl3 = rl_cTP(c,T3, P3);
s3 = s_crT(c, rl3, T3);
h3 = h_crT(c, rl3, T3);

% State 4: Air flashed to 1 bar
P4 = 1e5;
[T4, q4, V4, y4, x4, rg4, rf4] = Flash_zhP(c, h3, P4);
s4l = s_crT(x4, rf4, T4);
s4v = s_crT(y4, rg4, T4);
s4 = s4v*q4 + s4l*(1-q4);

%% Now iterate over output composition until column top matches air flash

% Column Parameters
n_trays = 3;
reboil_quality = 0.7;
P_col = 1e5;

% NR Parameters
imax = 100; i = 0;
x_inc = 1e-3;
err_tol = 1e-3;

% Make a guess
x_out(N2) = 0.9;
x_out(O2) = 1 - x_out(N2);

while true
    i = i + 1;
    
    [~, vapout(1), liqin(1)] = Reboiler_xqP(x_out,reboil_quality,P_col);
    for i = 2:n_trays+1
        [vapout(i), liqin(i)] = upTray(vapout(i-1), liqin(i-1));
    end
    err = abs( liqin(end).c(N2) - x4(N2));

    % Do NR
    % High
    x_high(N2) = x_out(N2) + x_inc;
    x_high(O2) = 1 - x_high(N2);
    [Q, vapout_high(1), liqin_high(1)] = Reboiler_xqP(x_high,reboil_quality,P_col);
    for i = 2:n_trays+1
        [vapout_high(i), liqin_high(i)] = upTray(vapout_high(i-1), liqin_high(i-1));
    end
    err_high = abs( liqin_high(end).c(N2) - x4(N2));
    % Low
    x_low(N2) = x_out(N2) - x_inc;
    x_low(O2) = 1 - x_low(N2);
    [Q, vapout_low(1), liqin_low(1)] = Reboiler_xqP(x_low,reboil_quality,P_col);
    for i = 2:n_trays+1
        [vapout_low(i), liqin_low(i)] = upTray(vapout_low(i-1), liqin_low(i-1));
    end
    err_low = abs( liqin_low(end).c(N2) - x4(N2));

    % Increment x_out
    dedx = (err_high - err_low)/(2*x_inc);
    x_out(N2) = x_out(N2) - err/dedx;
    x_out(O2) = 1 - x_out(N2);

    if i>imax
        disp("Failed to Converge Column")
        break
    elseif err < err_tol
        disp("Successfully solved for outlet composition")
        break
    end

end

x_out

% Assemble Data for Plotting
% Data from flashing air process
T_cycle = [T1, T1, T3, T4];
s_cycle = [s1, s2, s3, s4];
x_cycle = [c(N2), c(N2), c(N2), q4*y4(N2)+(1-q4)*x4(N2)]; % overall composition
T_tray = zeros(1,n_trays+1);
x_tray = zeros(1,n_trays+1);
s_tray = zeros(1,n_trays+1);
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


% Data for average
% Reboiler
% Tie Lines
y_tie(1) = vapout(1).c(N2);
x_tie(1) = x_out(N2);
T_tie(1) = vapout(1).T;
sv = s_crT(vapout(1).c, vapout(1).r, T_tie(1));
sv_tie(1) = sv;
r_xout = rl_cTP(x_out, T_tie(1), P_col);
sl = s_crT(x_out, r_xout, T_tie(1));
sl_tie(1) = sl;

% Process Lines
T_tray(1) = vapout(1).T;
x_tray(1) = reboil_quality*vapout(1).c(N2) + (1-reboil_quality)*x_out(N2);
s_tray(1) = sv*reboil_quality + sl*(1-reboil_quality);

for i = 2:length(vapout)
    % Tie Lines
    y_tie(i) = vapout(i).c(N2);
    x_tie(i) = liqin(i-1).c(N2);
    T_tie(i) = vapout(i).T;
    sv = s_crT(vapout(i).c, vapout(i).r, T_tie(i));
    sv_tie(i) = sv;
    sl = s_crT(liqin(i-1).c, liqin(i-1).r, T_tie(i));
    sl_tie(i) = sl;
    
    % Process Lines
    T_tray(i) = vapout(i).T;
    qual_tray = vapout(i).mdot/(vapout(i).mdot+liqin(i-1).mdot);
    x_tray(i) = qual_tray*vapout(1).c(N2) + (1-qual_tray)*liqin(i-1).c(N2);
    s_tray(i) = sv*qual_tray + sl*(1-qual_tray);
end

T_plot = [T_cycle flip(T_tray)];
x_plot = [x_cycle flip(x_tray)];
s_plot = [s_cycle flip(s_tray)];
% Start Plotting, include vapor dome data
% Load a file to save recalculating things up to here.

% We will build the surface as we go.  Start a figure and put it on hold.
figure(1)
clf
hold on
ylabel('Nitrogen Mole Fraction','rotation',0)
xlabel('Specific Entropy (kJ/kg-K)','rotation',0)
zlabel('Temperature (K)')
view([-10 40])
axis([2 7 0 1 60 300])
grid on
drawnow

% Plot Generated Data
plot3(s_plot/1e3,x_plot,T_plot,'-or','LineWidth',2);
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