% ME362B Problem Set 7 - Problem 2
% Winter 2023
% Andy Huynh

% Part 2ai
x = 0:0.005:0.15;
d = 0.01 + 2*x*tand(5);
A = pi*d.^2/4;

gamma = 5/3;

for i = 1:length(A)
    f = @(M) (((gamma+1)/2)^(-(gamma+1)/(2*(gamma-1)))*(1+(gamma-1)/2*M^2)^((gamma+1)/(2*(gamma-1)))/M-(A(i)/A(1)));
    Mach(i) = fzero(f,2);
end

% figure(1)
% hold on
% plot(x*100,Mach)
% xlabel('Distance Downstream(cm)')
% ylabel('Mach Number')
% xticks([0:1:15])
% hold off

% Part 2aii
Tt = 3000;
T = Tt*(1+((gamma-1)/2)*Mach.^2).^-1;

% figure(2)
% hold on
% plot(x*100,T)
% xlabel('Distance Downstream(cm)')
% ylabel('T_{tr} (K)')
% xticks([0:1:15])
% hold off

% Part 2aiii
Pt = 200; % atm
P = Pt*(T/Tt).^((gamma)/(gamma-1));

% figure(3)
% hold on
% plot(x*100,P)
% xlabel('Distance Downstream(cm)')
% ylabel('P (atm)')
% xticks([0:1:15])
% hold off

% Part 2aiv
R = 2076.9; % J/kg*K
a = (gamma*R*T).^0.5;
V = Mach.*a;

% figure(4)
% hold on
% plot(x*100,V)
% xlabel('Distance Downstream(cm)')
% ylabel('Velocity (m/s)')
% xticks([0:1:15])
% hold off


t = x./V;
% Teq = 3259;
%  
% x = 0:1/100:50/100;
% V1 = 5038; % velocity after shock
% t = x./V1;
%  
% T_tr(1) = 4047;
% T_v(1) = 140;
%  
% for i = 1:50
%     tau(i) = exp(220*(T_tr(i)^(-1/3)-0.015*14^0.25)-18.42)/P;
%     %f = @(Tb) (t(i)-t(i+1))/(log((1/(exp(3395/Teq)-1)-1/(exp(3395/Tb)-1))/(1/(exp(3395/Teq)-1)-1/(exp(3395/T_v(i))-1)))-tau(i);
%     f = @(Tb) (t(i)-t(i+1))/log((1/(exp(3395/Teq)-1) - 1/(exp(3395/Tb)-1)) / (1/(exp(3395/Teq)-1) - 1/(exp(3395/T_v(i))-1)))-tau(i);
%     T_v(i+1) = fzero(f,200);
%     f2 = @(T) (5/2)*(1.38e-23)*T+(1.38e-23)*3395/(exp(3395/T_v(i+1))-1)-1.396e-19;
%     T_tr(i+1) = fzero(f2,4000);
% end