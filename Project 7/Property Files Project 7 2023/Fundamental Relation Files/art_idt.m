function art = art_idt(i,d,t)
% Return the derivative of the residual normalized Helmholtz function of
% species i with respect to tau at a specified value of delta (d), and tau (t).

global FR_Npoly FR_Nexp FR_Ngaus FR_Nnonan
global FR_N FR_d FR_t FR_c
global FR_eta FR_beta FR_gamma FR_epsilon
global FR_a FR_b FR_NAbeta FR_A FR_B FR_C FR_D

% Build up the residual by the types of terms.
art = 0;
k = 1;

% Start with the polynomial terms:
for(j=1:1:FR_Npoly(i))
    art = art + FR_t(k,i)*FR_N(k,i)*d^FR_d(k,i)*t^(FR_t(k,i)-1);
    k = k+1;
end

% Add the exponential terms:
for(j=1:1:FR_Nexp(i))
    art = art + FR_t(k,i)*FR_N(k,i)*d^FR_d(k,i)*t^(FR_t(k,i)-1)*exp(-d^FR_c(k,i));
    k = k+1;
end

% Add the Gaussian terms:
for(j=1:1:FR_Ngaus(i))
    art = art + ...
        (...
        FR_t(k,i)*FR_N(k,i)*d^FR_d(k,i)*t^(FR_t(k,i)-1)...
        - 2*FR_beta(k,i)*(t-FR_gamma(k,i))*FR_N(k,i)*d^FR_d(k,i)*t^FR_t(k,i)...
        )...
        *exp(-FR_eta(k,i)*(d-FR_epsilon(k,i))^2-FR_beta(k,i)*(t-FR_gamma(k,i))^2);
    k = k+1;
end

% Add the nonanalytic terms:
for(j=1:1:FR_Nnonan(i))
    Theta = (1-t) + FR_A(k,i)*((d-1)^2)^(1/(2*FR_NAbeta(k,i)));
    Delta = Theta^2 + FR_B(k,i)*((d-1)^2)^FR_a(k,i);
    Psi   = exp(-FR_C(k,i)*(d-1)^2-FR_D(k,i)*(t-1)^2);

    % Don't compute the non-analytic terms if sitting on critical.
    % Their limit as you approach critical is zero.
    if((Theta == 0)&&(Delta == 0)&&(Psi == 1))
        art = art + 0;
    else
        tTheta = -1;
        tDelta = 2*Theta*tTheta;
        tPsi   = -2*FR_D(k,i)*(t-1)*exp(-FR_C(k,i)*(d-1)^2-FR_D(k,i)*(t-1)^2);

        art = art + FR_b(k,i)*FR_N(k,i)*Delta^(FR_b(k,i)-1)*tDelta*d*Psi...
            + FR_N(k,i)*Delta^FR_b(k,i)*d*tPsi;
    end
    k = k+1;
end
