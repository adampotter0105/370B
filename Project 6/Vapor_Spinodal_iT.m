function [vap_spin] = Vapor_Spinodal_iT(T)
%Vapor_Spinodal_iT returns the location/density (kg/m^3) of the vapor
%numerical spinodals (outermost locations at which dP/dr goes to zero) for
%any given temperature (K).

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

% Set limits and step size.
rmax = rcrit_i(ispecies);
rmin = rgtrip_i(ispecies)/2;

% Setting the hyperparameters.
m = 100; % number of samples
m_elite = 10; % number of elite samples
covar = 5*eye(1); % covariance matrix

% Starting point average.
i = 1;
r_list = linspace(rmin,rmax,m*10);
starting_dPdr = zeros(1,m);

for n = length(r_list):-1:1
    if d2Pdr2_irT(ispecies,r_list(n),T) > 0
        r_list(n) = []; % discard inner minimum
    end
end
for r = 1:length(r_list)
    starting_dPdr(i) = abs(dPdr_irT(ispecies,r_list(r),T));
    i = i+1;
end
[out,idx] = sort(starting_dPdr); % get index of sorted samples
elite_index = idx(1:m_elite); % retains the best m_elite samples' indices
elites = r_list(elite_index); % lists elite samples
avg = mean(elites); % calculates new mean of elite samples

dPdr_target = 0; % target minimum
dPdr_best = 9999; % initialize variable > tolerance

% Cross-entropy method.
% Samples that fall outside rmin or rmax are discarded (constraints).
while (dPdr_best-dPdr_target) > 0.01 % set tolerance
   samples =  mvnrnd(avg,covar,m); % creates multivariate normal distribution
   for k = m:-1:1
       if (samples(k) > rmax || samples(k) < rmin || d2Pdr2_irT(ispecies,samples(k),T) > 0)
           samples(k) = []; % discard infeasible samples
       end
   end
   for s = 1:length(samples)
       dPdr(s) = abs(dPdr_irT(ispecies,samples(s),T)); % evaluates dPdr with samples values
   end
   [out,idx] = sort(dPdr); % get index of sorted samples
   elite_index = idx(1:m_elite); % retains the best m_elite samples' indices
   elites = samples(elite_index); % lists elite samples
   avg = mean(elites); % calculates new mean of elite samples
   covar = cov(elites); % calculates new covariance matrix of elite samples
   dPdr_best = abs(dPdr_irT(ispecies,avg,T)); % returns dPdr using the last mean value
   r_best = min(elites); % modification to cross-entropy method to get absolute minimum
end

vap_spin = r_best; % best r thus far

end