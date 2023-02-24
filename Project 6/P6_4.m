% ME370B: Modeling and Advanced Concepts
% Project 6 - Part 4
% Dongwon Ka

% Make a pressure-specific volume isotherm diagram.  
% This file is useful for exploring how the P_vT function behaves in various 
% regions (triple line, critical point, etc.).
% C.F. Edwards, 2-11-12

% Provide access to support files via the Matlab path.
addpath 'Fundamental_Relation_Files' 
addpath 'Fundamental_Relation_Data'
addpath 'Setup_Files'
addpath 'Property_Files'

% Clean up and get ready to go.
clear all
format compact
fprintf('\n************************************************************\n')

% Set up the basic storage and load the FR files.
Setup_Props_i;

% Set which of the loaded species you want to work with.  You might want to
% change this around to see what the plots look like for other species.
ispecies = nH2;

%% Project 6 - Problem 4

% Pressure and enthalpy list for isobars
Plist = [0.1 0.2 0.5 1 2 5 10 20 50 100 200 500]*101325;
hlist = [0.35831 0.62672 0.89313 1.1605 1.428 1.6954 1.9628 2.2302 2.4976 2.765 3.0324 3.2998 3.5672 3.8346 4.1021 4.3695]*1e6;
%hlist = [1.1605 1.428 1.6954 1.9628 2.2302 2.4976 2.765 3.0324 3.2998 3.5672 3.8346 4.1021 4.3695]*1e6;

% temperature list
Tmin  = 13.96;
Tmax  = 310;
Tlist = linspace(Tmin,Tmax,100);

% triple line
num_Ttri  = 13.95;
num_rltri = 76.92;
num_rvtri = 0.145;
sltri     = s_irT(ispecies, num_rltri, num_Ttri); 
svtri     = s_irT(ispecies, num_rvtri, num_Ttri); 
T_tri     = num_Ttri*ones(100);
s_tri     = linspace(sltri,svtri,100);

fprintf("Generating Vapor Dome \n")
% vapor dome 
T_dome_x = linspace(1.05*Ttrip_i(ispecies), 0.99*Tcrit_i(ispecies), 15);    
s_dome_liq = NaN(1,length(T_dome_x));
s_dome_vap = NaN(1,length(T_dome_x));
for i = 1:length(T_dome_x)
    T = T_dome_x(i);
    [~, rl, rv] = Saturation_iT_lookup(T);
    s_dome_liq(i) = s_irT(ispecies, rl, T);
    s_dome_vap(i) = s_irT(ispecies, rv, T);
end

fprintf("Generating Isobars Under Dome \n")
% isobars under the dome (horizontal lines)
for j=1:1:length(Plist)
    P=Plist(j);
    if P < Pcrit_i(ispecies)
        [T, rl, rv] = Saturation_iP_lookup(P);
        T_dome(j)  = T;
        rv_list(j) = rv;
        s_dome_max = s_irT(ispecies, rl, T);
        s_dome_min = s_irT(ispecies, rv, T);
        s_dome(:,j)     = linspace(s_dome_min, s_dome_max, 100);
    end
end

fprintf("Generating Isobars in Vapor Phase \n")
% isobars vapor
for j=1:1:length(Plist)
    P=Plist(j)
    for i=1:1:length(Tlist)
        r = rv_iTP(ispecies, Tlist(i), P);
        s(i,j) = s_irT(ispecies, r, Tlist(i));
        if P < Pcrit_i(ispecies) %%%% To remove overlap the line under the dome
            if r > rv_list(j)
                s(i,j) =NaN;
            end
        end
    end
end

% isobars liquid (left part)
fprintf("Generating Isobars in Liquid Phase \n")
for j=1:1:length(Plist)
    P=Plist(j)
    for i=1:1:length(Tlist)
        if Tlist(i) < Tcrit_i(ispecies)
            T_left(i) = Tlist(i);
            r = rl_iTP(ispecies, Tlist(i), P);
            s_left(i,j) = s_irT(ispecies, r, Tlist(i));
        end
    end
end

fprintf("Generating Isenthalps \n")
% Isenthalps
Pmax = 500*101325;
Pmin = 0.1*101325;
Plist2 = linspace(Pmin, Pmax, 50);
T_for_hlines = NaN(length(Plist2),length(hlist));
s_for_hlines = NaN(length(Plist2),length(hlist));

for j=1:1:length(hlist)
    h=hlist(j)
    for i=1:1:length(Plist2)
        P = Plist2(i);
        [r, T, rf, rg]=rT_ihP(ispecies,h,P);
        T_for_hlines(i,j) = T;
        s_for_hlines(i,j) = s_irT(ispecies, r, T);

    end
end

fprintf("Starting Plotting \n")
%% Plotting
figure()
hold on
% plot for isobars (vapor)
for j=1:length(Plist)
    plot(s(:,j)/1e3,Tlist,'r')
end

% plot for isobars (liquid)
for j=1:length(Plist)
    plot(s_left(:,j)/1e3, T_left,'r')
end

% plot for isobar under the dome
for j=1:1:length(Plist)
    P=Plist(j);
    if P < Pcrit_i
        plot(s_dome(:,j)/1e3,T_dome(j)*ones(100), 'r')
    end
end

% plot for isenthalps w/o three low ethalpy
for j=1:length(hlist)
    plot(s_for_hlines(:,j)/1e3,T_for_hlines(:,j),'b')
end

% plot for triple line
plot(s_tri/1e3,T_tri,'k--','LineWidth',3)

%plot for vapor dome
plot(s_dome_liq/1e3, T_dome_x,'k-','LineWidth',3)
plot(s_dome_vap/1e3, T_dome_x,'k-','LineWidth',3)
hold off

ylim([0 350])
xlim([0 100])
plotfixer