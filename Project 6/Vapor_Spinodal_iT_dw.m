function [vap_spin] = Vapor_Spinodal_iT(ispecies, T)
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
ispecies = ispecies;

% Set limits and step size.
rmin = rgtrip_i(ispecies);
rmax = rcrit_i(ispecies); 

% Density range.
pts    = 100;
r_list = linspace(rmin,rmax,pts);

dPdr_list = 100000*ones(pts,1);
for i=1:pts
    if d2Pdr2_irT(ispecies, r_list(i),T)<0  % negative
        dPdr_list(i) = dPdr_irT(ispecies,r_list(i),T);
    end
end

for i=1:pts
    if abs(dPdr_list(i)) == min(abs(dPdr_list))
        vap_spin = r_list(i);
    end
end

