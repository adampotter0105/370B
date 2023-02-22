function [sat_pres, r_sat_liq, r_sat_vap] = Saturation_iT(T)
%Liquid_Spinodal_iT returns the location/density (kg/m^3) of the liquid
%numerical spinodals (outermost locations at which dP/dr goes to zero) for
%any given temperature (K).
%Dongwon Ka

% Provide access to support files via the Matlab path.
addpath 'Fundamental_Relation_Files'
addpath 'Fundamental_Relation_Data'
addpath 'Setup_Files'
addpath 'Property_Files'

% Set up the basic storage and load the FR files.
Setup_Props_i;

% Set which of the loaded species you want to work with.  You might want to
% change this around to see what the plots look like for other species.
ispecies = nH2;

% chemical potential min_spec_vol to liquid_spinodal point
% Set limits and step size.
rmax_1 = rupper_i(ispecies);
rmin_1 = Liquid_Spinodal_iT(T); 

% Density range.
pts      = 1000;
r_list_1 = linspace(rmin_1,rmax_1,pts);
for j=1:pts
    p_list_1(j)  = P_irT(ispecies, r_list_1(j), T);
    mu_list_1(j) = mu_irT(ispecies, r_list_1(j), T);
end

% chemical potential vapor_spinodal point to max_spec_vol
rmax_2   = Vapor_Spinodal_iT(T);
rmin_2   = rgtrip_i(ispecies);

% Density range
r_list_2 = linspace(rmin_2,rmax_2,pts);
for j=1:pts
    p_list_2(j)  = P_irT(ispecies, r_list_2(j), T);
    mu_list_2(j) = mu_irT(ispecies, r_list_2(j), T);
end

% Find the cross-point
test = 10000000;
num_r_list_1 = 0;
num_r_list_2 = 0;

for i=1:pts
    for j=1:pts
        if sqrt((p_list_2(i)-p_list_1(j))^2+(mu_list_2(i)-mu_list_1(j))^2) < test
            test = sqrt((p_list_2(i)-p_list_1(j))^2+(mu_list_2(i)-mu_list_1(j))^2);
            num_r_list_1 = j;
            num_r_list_2 = i;                
        end
    end
end

% output
sat_pres  = P_irT(ispecies, r_list_1(num_r_list_1),T);
r_sat_liq = r_list_1(num_r_list_1);
r_sat_vap = r_list_2(num_r_list_2);