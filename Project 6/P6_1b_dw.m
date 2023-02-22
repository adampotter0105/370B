% ME370B: Modeling and Advanced Concepts
% Project 6 - Part 1b
% Dongwon Ka

% Make a chemical potential-pressure isotherm diagram.  
% This file is useful for exploring how the mu_irT function behaves in various 
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
ispecies = nH2

% Put a few isotherms on the mu-P plot.
Tlist = [...
    Ttrip_i(ispecies)... % get ispecies-th element of Ttrip_i vector
    0.10*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.25*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.50*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.75*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.90*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.95*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    1.00*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    ]



% Preallocate storage...
% Prisotherm  = zeros(length(Tlist),steps+1);
% risotherm   = zeros(length(Tlist),steps+1);
% muisotherm  = zeros(length(Tlist),steps+1);

for j=1:1:length(Tlist)
    T = Tlist(j)
    % Set limits and step size.
    rmax_1 = rupper_i(ispecies);
    rmin_1 = Liquid_Spinodal_iT(ispecies,T); 
    steps = 1000;
    dr = (rmax_1-rmin_1)/steps;
    i = 1;
    for r=rmin_1:dr:rmax_1+dr
        r;
        Prisotherm_1(j,i) = P_irT(ispecies,r,T);
        muisotherm_1(j,i) = mu_irT(ispecies,r,T);
        risotherm_1(j,i)  = r;
        i = i+1;
    end

    rmax_2 = Vapor_Spinodal_iT(ispecies, T);
    rmin_2 = rgtrip_i(ispecies);
    dr = (rmax_2-rmin_2)/steps;
    i = 1;
    for r=rmin_2:dr:rmax_2+dr
        r;
        Prisotherm_2(j,i) = P_irT(ispecies,r,T);
        muisotherm_2(j,i) = mu_irT(ispecies,r,T);
        risotherm_2(j,i)  = r;
        i = i+1;
    end

end



mucrit = mu_irT(ispecies,rcrit_i(ispecies),Tcrit_i(ispecies)); % (J/kmol)

Tmin        = Ttrip_i(ispecies);
Tmax        = 0.999*Tcrit_i(ispecies);
T_steps     = 30;
T_vector_1  = linspace(Tmin,Tmax,T_steps);
%T_vector_2  = linspace(Tmax*5/10,Tmax,T_steps*20/30);
T_vector    = [T_vector_1];
for i=1:T_steps
    T_vector(i)
    [sat_pres(i), r_sat_liq(i), r_sat_vap(i)] = Saturation_iT(ispecies,T_vector(i));
    mu_sat_pres(i)                            = mu_irT(ispecies, r_sat_liq(i),T_vector(i));
end


figure(1)
clf
hold on
% Put the isotherms on.
for j=1:1:length(Tlist)
    plot(Prisotherm_1(j,:)/1e6,muisotherm_1(j,:)/1e6,'b')
    plot(Prisotherm_2(j,:)/1e6,muisotherm_2(j,:)/1e6,'b')
end
plot(sat_pres/1e6,mu_sat_pres/1e6,'k')
%plot(Prisotherm(1,:)/1e6,muisotherm(1,:)/1e6,'r') % red curves indicating boundaries
%plot(Prisotherm(length(Tlist),:)/1e6,muisotherm(length(Tlist),:)/1e6,'r')
%legend('Critical Point','Triple Line')
hold off
xlabel('Pressure (MPa)')
ylabel('Chemical Potential (MJ/kmol)')
% Add some simple temperature labels.
for i=2:1:length(Tlist)
    text(Pcrit_i(ispecies)/1e6,(i-2.0)*(mucrit/1e6)/length(Tlist),...
        num2str(Tlist(i),'%.2f'))
end
i = 1;
text(Pcrit_i(ispecies)/1e6,(i-2.0)*(mucrit/1e6)/length(Tlist),...
    ['T = ',num2str(Tlist(i)),' K'])
xlim([-Pcrit_i(ispecies)/4e6 1.75*Pcrit_i(ispecies)/1e6]); xticks([0 0.5 1 1.5 2]);
ylim([-0.8 0.2]); yticks([-0.8 -0.6 -0.4 -0.2 0 0.2]);

% Gussy up the plot a little.
plotfixer