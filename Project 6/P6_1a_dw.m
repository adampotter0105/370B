% ME370B: Modeling and Advanced Concepts
% Project 6 - Part 1a
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
ispecies = nH2

% Put a few isotherms on the P-v plot.
Tlist = [...
    Ttrip_i(ispecies)...
    0.1*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.25*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.5*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.75*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.9*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.95*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    1.00*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    1.05*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    1.25*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    1.50*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    2.0*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.1*(Tupper_i(ispecies)-Tcrit_i(ispecies)) + Tcrit_i(ispecies) ...
    0.4*(Tupper_i(ispecies)-Tcrit_i(ispecies)) + Tcrit_i(ispecies) ...
    Tupper_i(ispecies)...
    ]

% Set limits and step size.
rmax = rupper_i(ispecies);
rmin = rgtrip_i(ispecies)/2;
steps = 2000;
dr = (rmax-rmin)/steps;
% Preallocate storage...
Prisotherm  = zeros(length(Tlist),steps+1);
risotherm   = zeros(length(Tlist),steps+1);

%%
for j=1:1:length(Tlist)
    T = Tlist(j)
    i = 1;
    for r=rmin:dr:rmax+dr
        r;
        Prisotherm(j,i) = P_irT(ispecies,r,T);
        risotherm(j,i) = r;
        i = i+1;
    end

end

for j=1:1:length(Tlist)
    T = Tlist(j);
    if T<Tcrit_i
        r_LS_on_isotherm(j) = Liquid_Spinodal_iT(ispecies,T);
        p_LS_on_isotherm(j) = P_irT(ispecies,r_LS_on_isotherm(j),T);
        r_VS_on_isotherm(j) = Vapor_Spinodal_iT(ispecies,T);
        p_VS_on_isotherm(j) = P_irT(ispecies,r_VS_on_isotherm(j),T);
    
        i=1;
        while (1)
            if risotherm(j,i) > r_VS_on_isotherm(j)
                min_r_under_dome = risotherm(j,i);
                num_min_r        = i;
                break;
            end
            i=i+1;
        end    
        i=steps;
        while (1)
            if risotherm(j,i) < r_LS_on_isotherm(j)
                max_r_under_dome = risotherm(j,i);
                num_max_r        = i;
                break;
            end
            i=i-1;
        end

        if j==1
            risotherm_1_1(:)  = risotherm(j,1:num_min_r-1);   risotherm_1_2(:) = risotherm(j,num_max_r+1:steps+2);
            Prisotherm_1_1(:) = Prisotherm(j,1:num_min_r-1); Prisotherm_1_2(:) = Prisotherm(j,num_max_r+1:steps+2);
        elseif j==2
            risotherm_2_1(:)  = risotherm(j,1:num_min_r-1);   risotherm_2_2(:) = risotherm(j,num_max_r+1:steps+2);
            Prisotherm_2_1(:) = Prisotherm(j,1:num_min_r-1); Prisotherm_2_2(:) = Prisotherm(j,num_max_r+1:steps+2);
        elseif j==3
            risotherm_3_1(:)  = risotherm(j,1:num_min_r-1);   risotherm_3_2(:) = risotherm(j,num_max_r+1:steps+2);
            Prisotherm_3_1(:) = Prisotherm(j,1:num_min_r-1); Prisotherm_3_2(:) = Prisotherm(j,num_max_r+1:steps+2);
        elseif j==4
            risotherm_4_1(:)  = risotherm(j,1:num_min_r-1);   risotherm_4_2(:) = risotherm(j,num_max_r+1:steps+2);
            Prisotherm_4_1(:) = Prisotherm(j,1:num_min_r-1); Prisotherm_4_2(:) = Prisotherm(j,num_max_r+1:steps+2);
        elseif j==5
            risotherm_5_1(:)  = risotherm(j,1:num_min_r-1);   risotherm_5_2(:) = risotherm(j,num_max_r+1:steps+2);
            Prisotherm_5_1(:) = Prisotherm(j,1:num_min_r-1); Prisotherm_5_2(:) = Prisotherm(j,num_max_r+1:steps+2);
        elseif j==6
            risotherm_6_1(:)  = risotherm(j,1:num_min_r-1);   risotherm_6_2(:) = risotherm(j,num_max_r+1:steps+2);
            Prisotherm_6_1(:) = Prisotherm(j,1:num_min_r-1); Prisotherm_6_2(:) = Prisotherm(j,num_max_r+1:steps+2);      
        elseif j==7
            risotherm_7_1(:)  = risotherm(j,1:num_min_r-1);   risotherm_7_2(:) = risotherm(j,num_max_r+1:steps);
            Prisotherm_7_1(:) = Prisotherm(j,1:num_min_r-1); Prisotherm_7_2(:) = Prisotherm(j,num_max_r+1:steps);          
        end
    end
end

Tmin        = Ttrip_i(ispecies);
Tmax        = 0.999*Tcrit_i(ispecies);  % We should find the numerical critical point from 3(a)
T_steps     = 30;
T_vector_1  = linspace(Tmin,Tmin+(Tmax-Tmin)*0.5/10,T_steps*10/30);
T_vector_2  = linspace(Tmin+(Tmax-Tmin)*0.5/10,Tmin+(Tmax-Tmin)*9/10,T_steps*5/30);
T_vector_3  = linspace(Tmin+(Tmax-Tmin)*9/10,Tmax,T_steps*15/30);
T_vector    = [T_vector_1 T_vector_2 T_vector_3];
for i=1:T_steps
    T_vector(i)
    r_LS(i) = Liquid_Spinodal_iT(ispecies,T_vector(i));
    p_LS(i) = P_irT(ispecies,r_LS(i),T_vector(i));
    r_VS(i) = Vapor_Spinodal_iT(ispecies,T_vector(i));
    p_VS(i) = P_irT(ispecies,r_VS(i),T_vector(i));
    [Sat_p(i), r_sat_liq(i), r_sat_vap(i)] = Saturation_iT(ispecies, T_vector(i));
end

figure(1)
clf
% Put the critical point and ends of the triple line on.
semilogx(1/rcrit_i(ispecies),Pcrit_i(ispecies)/1e6,'kd')
hold on
semilogx([1/rftrip_i(ispecies) 1/rgtrip_i(ispecies)],...
   [Ptrip_i(ispecies) Ptrip_i(ispecies)]/1e6,'ko-')
% Put the isotherms on.
semilogx(1./risotherm_1_1,Prisotherm_1_1/1e6,'r')
semilogx(1./risotherm_2_1,Prisotherm_2_1/1e6,'b')
semilogx(1./risotherm_3_1,Prisotherm_3_1/1e6,'b')
semilogx(1./risotherm_4_1,Prisotherm_4_1/1e6,'b')
semilogx(1./risotherm_5_1,Prisotherm_5_1/1e6,'b')
semilogx(1./risotherm_6_1,Prisotherm_6_1/1e6,'b')
semilogx(1./risotherm_7_1,Prisotherm_7_1/1e6,'b')
semilogx(1./risotherm_1_2,Prisotherm_1_2/1e6,'r')
semilogx(1./risotherm_2_2,Prisotherm_2_2/1e6,'b')
semilogx(1./risotherm_3_2,Prisotherm_3_2/1e6,'b')
semilogx(1./risotherm_4_2,Prisotherm_4_2/1e6,'b')
semilogx(1./risotherm_5_2,Prisotherm_5_2/1e6,'b')
semilogx(1./risotherm_6_2,Prisotherm_6_2/1e6,'b')
semilogx(1./risotherm_7_2,Prisotherm_7_2/1e6,'b')
for j=8:1:length(Tlist)
    semilogx(1./risotherm(j,:),Prisotherm(j,:)/1e6,'b')
end
% red curves indicating boundaries
%semilogx(1./risotherm(1,:),Prisotherm(1,:)/1e6,'r') 
semilogx(1./risotherm(length(Tlist),:),Prisotherm(length(Tlist),:)/1e6,'r')
% spinodal
semilogx(1./r_LS,p_LS/1e6,'k--')
semilogx(1./r_VS,p_VS/1e6,'k--')
% vapor dome
semilogx(1./r_sat_liq,Sat_p/1e6,'k')
semilogx(1./r_sat_vap,Sat_p/1e6,'k')
%

legend('Critical Point','Triple Line')
hold off
xlabel('Specific Volume (m^3/kg)')
ylabel('Pressure (MPa)')
% Add some simple temperature labels.
for i=2:1:length(Tlist)-1
    text(2*1/rcrit_i(ispecies),3*(i-5.5)*(Pcrit_i(ispecies)/1e6)/length(Tlist),...
        num2str(Tlist(i),'%.2f'))
end
i = 1;
text(2*1/rcrit_i(ispecies),3*(i-5.5)*(Pcrit_i(ispecies)/1e6)/length(Tlist),...
    num2str(Tlist(i)),'Color','r')
i = length(Tlist);
text(2*1/rcrit_i(ispecies),3*(i-5.5)*(Pcrit_i(ispecies)/1e6)/length(Tlist),...
    ['T = ',num2str(Tlist(i)),' K'],'Color','r')
%ylim([-Pcrit_i(ispecies)/1e6 2*Pcrit_i(ispecies)/1e6])
ylim([-0.15 2])
xlim([0.01 10])
% Gussy up the plot a little.
plotfixer