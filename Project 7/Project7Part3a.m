% Project 7 Part 3a

addpath 'Fundamental Relation Files'
addpath 'Fundamental Relation Data'
addpath 'Mixture Models'
addpath 'Setup Files' 
addpath 'Property Files'
addpath 'Procedure Files'

clear all
format compact
fprintf('\n************************************************************\n')

% Load a file to save recalculating things up to here.
load Tsx_Data

N = 2;
Setup_Air_Props;


% Molar masses:
MW_N2 = 28.0134;    % (g/mol)
MW_O2 = 31.999;     % (g/mol)

% Set an air composition.
c = zeros(N,1);
% Use binary engineering air.
c(O2) = 0.21;
c(N2) = 1 - c(O2);

% Get the inflection point for this composition.
[Tinfl, rinfl]   = Pr_Inflection_c(c);
Pinfl           = P_crT(c,rinfl,Tinfl);

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
% Tell the user the composition planes to be shown.
xN2_planes = clist

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
for j=1:1:1
    plot3(sdewSpline(j,:)/1e3,xN2Spline,TdewSpline(j,:),'b-','LineWidth',2)
    plot3(sbubSpline(j,:)/1e3,xN2Spline,TbubSpline(j,:),'b-','LineWidth',2)
    drawnow
end

% Tell the user the pressure surfaces to be shown.
P_surfaces = Plist;

% Add isobars, just for P = 1
for i=1:1:NCS
    % Do the points below the dome.
    for j=1:1:1
        plot3(sisoPlow(:,j,i)/1e3,xN2(i)*ones(length(sdome)/2,1),TisoPlow(:,j,i),'b-','LineWidth',2)
        plot3(sisoPhigh(:,j,i)/1e3,xN2(i)*ones(length(sdome)/2,1),TisoPhigh(:,j,i),'b-','LineWidth',2)
        drawnow
    end
    % Now do the pressures that do not have dew points (above dome).
    for j=j+1:1:1
        plot3(sisoPlow(:,j,i)/1e3,xN2(i)*ones(length(sdome)/2,1),TisoPlow(:,j,i),'b-','LineWidth',2)
        plot3(sisoPhigh(:,j,i)/1e3,xN2(i)*ones(length(sdome)/2,1),TisoPhigh(:,j,i),'b-','LineWidth',2)
        drawnow
    end
end

Trange = linspace(60,300,10);
P_iso = 100000; % 1 bar
qual_mass = zeros(1,length(Trange));
V_mole = zeros(1,length(Trange));
N2_y_vap = zeros(1,length(Trange));
N2_x_liq = zeros(1,length(Trange));
rg_gas = zeros(1,length(Trange));
rf_liq = zeros(1,length(Trange));
s_vap = zeros(1,length(Trange));
s_liq = zeros(1,length(Trange));

for k = length(Trange):-1:1
    T = Trange(k)
    [qual_mass(k), V_mole(k), y_vap_hold, x_liq_hold, rg_gas(k), rf_liq(k)] = Flash_zTP(c,T,P_iso);
    N2_y_vap(k) = y_vap_hold(N2);
    N2_x_liq(k) = x_liq_hold(N2);
    s_vap(k) = s_crT(y_vap_hold, rg_gas(k),T);
    s_liq(k) = s_crT(x_liq_hold, rf_liq(k),T);
    plot3(s_vap(k)/1e3,N2_y_vap(k),Trange(k),'ro','LineWidth',2);
    plot3(s_liq(k)/1e3,N2_x_liq(k),Trange(k),'go','LineWidth',2);
    drawnow;
end

% plot3(s_vap/1e3,N2_y_vap,Trange,'ro','LineWidth',2);
% plot3(s_liq/1e3,N2_x_liq,Trange,'go','LineWidth',2);
% drawnow
hold off;
