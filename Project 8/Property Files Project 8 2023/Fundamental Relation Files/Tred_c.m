function T = Tred_c(c)
% Find the reducing T for this composition.

global Tcrit_i
global Zeta Beta Phi

T = 0;
for(i=1:1:length(c))
    T = T + c(i)*Tcrit_i(i);
end
for(i=1:1:length(c)-1)
    for(j=i+1:1:length(c))
        T = T + c(i)^Beta(i,j)*c(j)^Phi(i,j)*Zeta(i,j);
    end
end