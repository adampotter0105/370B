% ME370B: Modeling and Advanced Concepts
% Project 6 - Part 1b
% Andy Huynh

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

% Set limits and step size.
rmax = rupper_i(ispecies);
rmin = rgtrip_i(ispecies)/2;
steps = 2000;
dr = (rmax-rmin)/steps;
% Preallocate storage...
Prisotherm  = zeros(length(Tlist),steps+1);
risotherm   = zeros(length(Tlist),steps+1);
muisotherm  = zeros(length(Tlist),steps+1);
for j=1:1:length(Tlist)
    T = Tlist(j)
    i = 1;
    for r=rmin:dr:rmax+dr
        r;
        Prisotherm(j,i) = P_irT(ispecies,r,T);
        muisotherm(j,i) = mu_irT(ispecies,r,T);
        risotherm(j,i)  = r;
        i = i+1;
    end
end

mucrit = mu_irT(ispecies,rcrit_i(ispecies),Tcrit_i(ispecies)); % (J/kmol)

figure(1)
clf
hold on
% Put the isotherms on.
for j=1:1:length(Tlist)
    plot(Prisotherm(j,:)/1e6,muisotherm(j,:)/1e6,'b')
end
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
% plotfixer