function aedddd = aedddd_cdt(c,d,t)
% Return the 4th derivative of the excess contribution wrt delta
% at a specified value of composition (c), delta (d), and tau (t).

global Fmix Mix_i Mix_j Mix_N NMix

% Use the expression for the fit.

f = 0;
for(k=1:1:NMix)
    f = f + Mix_N(k)*Mix_i(k)*(Mix_i(k)-1)*(Mix_i(k)-2)*(Mix_i(k)-3)...
        *d^(Mix_i(k)-4)*t^Mix_j(k);
end

aedddd = 0;
for(i=1:1:length(c)-1)
    for(j=i+1:1:length(c))
        aedddd = aedddd + c(i)*c(j)*Fmix(i,j)*f;
    end
end