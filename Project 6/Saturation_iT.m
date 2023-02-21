function [P_sat,r_liq,r_vap] = Saturation_iT(T)
%Saturation_iT returns the saturation pressure (Pa), saturated liquid
%density, and saturated vapor density (kg/m^3) for given T (K) between the
%triple temperature and near-critical temperature.
    % Uses Maxwell's Equal-Area Rule

% Provide access to support files via the Matlab path.
addpath 'Fundamental_Relation_Files'
addpath 'Fundamental_Relation_Data'
addpath 'Setup_Files'
addpath 'Property_Files'

% Set up the basic storage and load the FR files.
Setup_Props_i;

% Set which of the loaded species you want to work with.  You might want to
% change this around to see what the plots look like for other species.
ispecies = nH2

% Setting the hyperparameters.
m = 100; % number of samples
m_elite = 10; % number of elite samples
covar = eye(1); % covariance matrix

% Starting point average.
liq_spin = Liquid_Spinodal_iT(T);
vap_spin = Vapor_Spinodal_iT(T);
P1 = P_irT(ispecies,liq_spin,T)/1e6;
P2 = P_irT(ispecies,vap_spin,T)/1e6;
avg = (P1+P2)/2; % (MPa)

dP_target = 0; % target area difference above/below P
dP_best = 9999; % initialize variable > tolerance

% Cross-entropy method.
while abs((dP_best-dP_target)) > 0.1 % set tolerance
   samples =  mvnrnd(avg,covar,m); % creates multivariate normal distribution
   for s = 1:length(samples)
       vmin = 1/rl_iTP(ispecies,T,samples(s)); %something is wrong here - I want to get the saturation v at specified P
       vmax = 1/rv_iTP(ispecies,T,samples(s));
       x = linspace(vmin,vmax,200);
       p = 1;
       for j = 1:length(x)
           Pvisotherm(p) = P_irT(ispecies,1/x(j),T)/1e6;
           p = p+1;
       end
       y = Pvisotherm - ones(1,200)*samples(s);
       dP(s) = abs(trapz(x,y));
   end
   [out,idx] = sort(dP); % get index of sorted samples
   elite_index = idx(1:m_elite); % retains the best m_elite samples' indices
   elites = samples(elite_index); % lists elite samples
   avg = mean(elites); % calculates new mean of elite samples
   covar = cov(elites); % calculates new covariance matrix of elite samples
   P_best = min(elites); % modification to cross-entropy method to get absolute minimum
   vmin = 1/rl_iTP(ispecies,T,avg);
   vmax = 1/rv_iTP(ispecies,T,avg);
   x = linspace(vmin,vmax,200);
   y = ones(1,200)*samples(s);
   dP_best = trapz(x,y); % returns dP using the last mean value
end

P_sat = P_best; % best P thus far (MPa)
r_liq = rl_iTP(ispecies,T,P_sat); % liquid density (kg/m^3) at P_sat, T
r_vap = rv_iTP(ispecies,T,P_sat); % vapor density (kg/m^3) at P_sat, T

end