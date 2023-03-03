% Adam Potter 3/2/23

addpath 'Fundamental Relation Files'
addpath 'Fundamental Relation Data'
addpath 'Mixture Models'
addpath 'Setup Files' 
addpath 'Property Files'
addpath 'Procedure Files'

clear all
format compact
fprintf('\n************************************************************\n')

%% Set up Air properties
N = 2;
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

% Set an air composition.
c(O2) = 0.21;
c(N2) = 1 - c(O2);

% Goal: T-s-x diagram: Start with air at ambient conditions, compress it 
% isothermally to 150 bar and then cool it at constant pressure to 150 K. 
% Flash the fluid (via throttle) from this state into a receiver (tank) at 
% 5 bar. Return the vapor from that receiver to ambient temperature by 
% isobaric heat transfer, and flash (throttle) the liquid to another 
% receiver maintained at a pressure of 1 bar. (The liquid that results from
% the second flash is the product and is sometimes referred to as 
% oxygen-enriched liquid air.) Return the vapor from that flash to ambient 
% temperature by isobaric heat transfer.

% State 1: Ambient air
% State 2: 150 bar, T_amb
% State 3: 150 bar, 150k (all vapor up till now)
% State 4: isenthalpic flash down to 5 bar 
% State 5: just vapor returns to 5 bar, T_amb
% State 6: isenthalpic flash down to 1 bar
% State 7: just vapor return to T_amb, 1bar

T_amb = 298;
P_amb = oneatm;

% State 1
rv = rv_cTP(c,T_amb, P_amb);
h1 = h_crT(c, rv, T_amb);
s1 = s_crT(c, rv, T_amb);
x1 = c;
T1 = T_amb;

% State 2
rv = rv_cTP(c,T_amb, 150*1e5);
h2 = h_crT(c, rv, T_amb);
s2 = s_crT(c, rv, T_amb);
x2 = c;
T2 = T_amb;

% State 3
rv = rv_cTP(c, 150, 150*1e5);
h3 = h_crT(c, rv, 150);
s3 = s_crT(c, rv, 150);
x3 = c;
T3 = 150;

% State 4
[~, ~, y, x, rg, rf, T] = Flash_zhP(c, h3, 5*1e5);
h4_l = h_crT(x, rf, T);
h4_v = h_crT(y, rg, T);
s4_l = s_crT(x, rf, T);
s4_v = s_crT(y, rg, T);
x4_l = x;
x4_v = y;
T4 = T;

% State 5
rv = rv_cTP(x, T_amb, 5*1e5);
h5 = h_crT(x, rv, T_amb);
s5 = s_crT(x, rv, T_amb);
x5 = x;
T5 = T_amb;

% State 6
[~, ~, y, x, rg, rf, T] = Flash_zhP(x, h5, 1e5);
h6_l = h_crT(x, rf, T);
h6_v = h_crT(y, rg, T);
s6_l = s_crT(x, rf, T);
s6_v = s_crT(y, rg, T);
x6_l = x;
x6_v = y;
T6 = T;

% State 7 (return to ambient vapor)
rv = rv_cTP(y, T_amb, P_amb);
h7 = h_crT(y, rv, T_amb);
s7 = s_crT(y, rv, T_amb);
x7 = y;
T7 = T_amb;


% Plot Process in T-s-x
figure(1)
hold on
xlabel('Specific Entropy (kJ/kg-K)','rotation',0)
ylabel('Nitrogen Mole Fraction','rotation',0)
zlabel('Temperature (K)')

% Collect Variables into arrays
x_vap = [x1(N2), x2(N2), x3(N2), x4_v(N2), x5(N2), x6_v(N2), x7(N2)];
s_vap = [s1, s2, s3, s4_v, s5, s6_v, s7];
T_vap = [T1, T2, T3, T4, T5, T6, T7];
x_liq1 = [x3(N2), x4_l(N2)];
s_liq1 = [s3, s4_l];
T_liq1 = [T3, T4];
x_liq2 = [x5(N1), x6_l(N2)];
s_liq2 = [s5, s6_l];
T_liq2 = [T5, T6];

%% Start Plotting, include vapor dome data
% Load a file to save recalculating things up to here.
load Tsx_Data

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

% Plot Cycle
plot3(s_vap/1e3, x_vap, T_vap, "ro--")
plot3(s_liq1/1e3, x_liq1, T_liq1, "ro--")
plot3(s_liq2, x_liq2, T_liq2, "r--")

% Tell the user the composition planes to be shown.
xN2_planes = clist;

% Show the P-rho critical locus on the plot.
plot3(sc/1e3,xN2,Tc,'kx','LineWidth',2)

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
for j=[1 3]
    plot3(sdewSpline(j,:)/1e3,xN2Spline,TdewSpline(j,:),'b-','LineWidth',2)
    plot3(sbubSpline(j,:)/1e3,xN2Spline,TbubSpline(j,:),'b-','LineWidth',2)
    drawnow
end

% Tell the user the pressure surfaces to be shown.
P_surfaces = Plist

% Add isobars.
for i=1:1:NCS
    % Do the points below the dome.
    for j=[1 3]
        plot3(sisoPlow(:,j,i)/1e3,xN2(i)*ones(length(sdome)/2,1),TisoPlow(:,j,i),'b-','LineWidth',2)
        plot3(sisoPhigh(:,j,i)/1e3,xN2(i)*ones(length(sdome)/2,1),TisoPhigh(:,j,i),'b-','LineWidth',2)
        drawnow
    end
end

% Finished with the dynamic plot.  Close it out.
hold off
improvePlot