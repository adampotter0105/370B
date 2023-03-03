% Adam Potter 3/2/23

addpath 'Fundamental Relation Files'
addpath 'Fundamental Relation Data'
addpath 'Mixture Models'
addpath 'Setup Files' 
addpath 'Property Files'
addpath 'Procedure Files'

%clear all
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


% State 1: Ambient air
% State 2: 1 bar, 60k

Tlist = linspace(60,300,400); % Need this many points to have 4 within the vapor dome
P = 1e5;

ylist = NaN(2,length(Tlist));
xlist = NaN(2,length(Tlist));
rglist = NaN(2,length(Tlist));
rflist = NaN(2,length(Tlist));
sl_list = NaN(1,length(Tlist));
sv_list = NaN(1,length(Tlist));
slist = NaN(1,length(Tlist));
c_list = repmat(c,length(Tlist),1);

% T_din = Dew_cP(c,P);
% T_dout = Bubble_cP(c,P);
% Commented out to run faster
T_din = 81.5557;
T_dout = 78.7277;

for i = 1:length(Tlist)
    T = Tlist(i)
    if T<=T_din && T>=T_dout
        [~, ~, ylist(:,i), xlist(:,i), rglist(:,i), rflist(:,i)] = Flash_zTP(c, T, P);
        sl_list(i) = s_crT(xlist(1,i), rflist(1,i), T);
        sv_list(i) = s_crT(ylist(1,i), rglist(1,i), T);
    elseif T>T_din % Vapor
        rv = rv_cTP(c, T, P);
        slist(i) = s_crT(c, rv, T);
    else % Liquid
        rl = rl_cTP(c, T, P);
        slist(i) = s_crT(c, rl, T);
    end
end


% Plot Process in T-s-x
figure(1)
hold on
xlabel('Specific Entropy (kJ/kg-K)','rotation',0)
ylabel('Nitrogen Mole Fraction','rotation',0)
zlabel('Temperature (K)')

%% Start Plotting, include vapor dome data
% Load a file to save recalculating things up to here.

% We will build the surface as we go.  Start a figure and put it on hold.
figure(1)
clf
hold on
ylabel('Nitrogen Mole Fraction','rotation',0)
xlabel('Specific Entropy (J/kg-K)','rotation',0)
zlabel('Temperature (K)')
view([-10 40])
axis([2 7 0 1 60 300])
grid on
drawnow

% Plot Tie Lines
n_ties = 4;
sl_list_cpy = sl_list;
for k = 1:n_ties
    [~, idx] = max(sl_list_cpy);
    plot3([sl_list(1,idx) sv_list(1,idx)]/1e3, [xlist(1,idx) ylist(1,idx)], [Tlist(idx) Tlist(idx)], "ro--")
    sl_list_cpy(idx) = NaN;
end

% Plot Cycle
plot3(slist/1e3, c_list(:,1), Tlist, "r")
plot3([slist(32) slist(37)]/1e3, c_list(1:2,1), [Tlist(32) Tlist(37)], "r")
plot3(sl_list(1,:)/1e3, xlist(1,:), Tlist, "r")
plot3(sv_list(1,:)/1e3, ylist(1,:), Tlist, "r")

load Tsx_Data
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
for j=1
    plot3(sdewSpline(j,:)/1e3,xN2Spline,TdewSpline(j,:),'b-','LineWidth',2)
    plot3(sbubSpline(j,:)/1e3,xN2Spline,TbubSpline(j,:),'b-','LineWidth',2)
    drawnow
end

% Tell the user the pressure surfaces to be shown.
P_surfaces = Plist

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
plotfixer