function a0t = a0t_idt(i,d,t)
% Return the derivative of the ideal gas contribution to the normalized 
% Helmholtz function with respect to tau at a specified value of 
% delta (d), and tau (t).

global O2
global FR_Npoly0 FR_Neinst FR_Nspec0
global FR_N0 FR_t0 FR_gamma0 FR_spec0

% Build up the residual by the types of terms.
a0t = 0;
k = 1;

% Start with the polynomial terms:
for(j=1:1:FR_Npoly0(i))
    if(FR_t0(k,i) == 0)
        a0t = a0t;
    else
        a0t = a0t + FR_t0(k,i)*FR_N0(k,i)*t^(FR_t0(k,i)-1);
    end
    k = k+1;
end

% Add the logarithmic terms:
a0t = a0t + FR_N0(k,i)/t;
k = k+1;
a0t = a0t + FR_N0(k,i)*log(t) + FR_N0(k,i);
k = k+1;

% Add the Einstein terms:
for(j=1:1:FR_Neinst(i))
    a0t = a0t + FR_N0(k,i)*FR_gamma0(k,i)*((1-exp(-FR_gamma0(k,i)*t))^(-1) - 1);
    k = k+1;
end

% If there are no special terms, return the answer now.
if(FR_Nspec0(i) == 0)
    return
end

% Add the special terms:
switch(i)
    case(O2)
        a0t = a0t + FR_N0(k,i)*(exp(FR_spec0(k,i)*t)-1)^(-1)...
            *exp(FR_spec0(k,i)*t)*FR_spec0(k,i);
        k = k+1;
        a0t = a0t + FR_N0(k,i)*(1+(2/3)*exp(-FR_spec0(k,i)*t))^(-1)...
            *(2/3)*exp(-FR_spec0(k,i)*t)*(-FR_spec0(k,i));
        return
    otherwise
        disp('Mistaken attempt to execute special terms in a0_idt')
        return
end
