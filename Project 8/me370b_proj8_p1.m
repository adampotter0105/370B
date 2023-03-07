clear all

%% Set up Air properties
N = 3;
Setup_Air_Props
Mair = 28.9586;                 % kg/kmol
Tmaxcondentherm = 132.6312;     % K
Pmaxcondentherm = 3.78502e6;    % Pa
rmaxcondentherm = 10.4477*Mair; % kg/m3
Tmaxcondenbar   = 132.6035;     % K
Pmaxcondenbar   = 3.7891e6;     % Pa
rmaxcondenbar   = 11.0948*Mair; % kg/m3
Tcritair        = 132.5306;     % K
Pcritair        = 3.7860e6;     % Pa
rcritair        = 11.8308*Mair; % kg/m3
% Bottom of dome conditions for air:
Tsolidair       = 59.75;        % K
Psolidair       = 5265;         % Pa
rsolidair       = 33.067*Mair;  % kg/m3
% Upper limit conditions for the model:
Tupper          = 870;          % K
Pupper          = 100e6;        % Pa
rupper          = rsolidair;
% Lower limit to stay just above solid air:
Tlower          = 60;           % K
% Molar masses
MW_O2 = 31.999; % (g/mol)
MW_N2 = 28.0134; % (g/mol)
MW_Ar = 39.948; % (g/mol)

Setup_Air_Props

% Functions for saturated liquids and vapors
% [T rl rv y] = Fast_Bubble_cP(x,P,varargin)
% [T rv rl x] = Fast_Dew_cP(y,P,varargin)
%         Tstart  = varargin{1};
%         xstart  = varargin{2};
%         rlstart = varargin{3};
%         rvstart = varargin{4};
% Functions for multi-component flash processes under the vapor dome
% [T x y rf rg] = Flash_zqP(z,q,P,varargin)
% [T q V y x rg rf] = Flash_zhP(z,h,P,varargin)
% [q V y x rg rf] = Flash_zTP(z,T,P,varargin)
%         T  = varargin{1};
%         y  = varargin{2};
%         rf = varargin{3};
%         rg = varargin{4};


%% Script for forward problem 
% % Inputs: Inlet composition, 
% x_in(N2) = 0.78; x_in(O2) = 0.21; x_in(Ar) = 0.01;
% qual_out = 0.8; % output mass quality vapor
% 
% % Reboiler Column Specifications
% ndot_out = 1; % kmol/s
% P_boil = 1e5; % Pressure in column known
% 
% % Interate on x_boil
% 
% x_boil = x_in;
% 
% % Find composition of streams using input composition
% [T_in, rl_in, ~, ~] = Fast_Bubble_cP(x_in,P_boil);
% [T_boil, rl_boil, rv_boil, y_boil] = Fast_Bubble_cP(x_boil,P_boil);
% h_in = h_crT(x_in, rl_in, T_in);
% h_vap = h_crT(y_boil, rv_boil, T_boil);
% h_out = h_crT(x_boil, rl_boil, T_boil);
% 
% % Find mass flow using quality of output
% mdot_out = ndot_out * sum([MW_N2, MW_O2, MW_Ar].*x_boil); % kmol*(kg/kmol/s) = kg/s
% mdot_vap  = mdot_out*(1/qual_out  - 1); % kg/s
% ndot_vap = mdot_vap / sum([MW_N2, MW_O2, MW_Ar].*y_boil); % kg/(kg/kmol/s) = kmol/s
% ndot_in = ndot_vap + ndot_out;
% 
% % First law balance
% Q = h_out*ndot_out + h_vap*ndot_vap - h_in*ndot_in;
% 
% % Output Q, [T_boil, rl_boil, x_boil], T_in

%% Call Function to Solve backward Problem

x_boil(O2) = 0.95; x_boil(N2) = 0.025; x_boil(Ar) = 0.025;
P_boil = 1e5;
qualities = 0:0.1:1;

x_in = NaN(3, length(qualities));
T_boil = NaN(1, length(qualities));
T_in = NaN(1, length(qualities));
Q = NaN(1, length(qualities));
for i = 1:length(qualities)
    [Q(i), T_boil(i), ~, ~, x_in(:,i), T_in(i)] = Reboiler_cqP(x_boil, qualities(i), P_boil);
end

% Plot Input Composition vs quality
figure(1)
plot(qualities, x_in(N2,:), "r-O")
hold on
plot(qualities, x_in(O2,:), "g-O")
plot(qualities, 10*x_in(Ar,:), "b-O")
plot([qualities(1) qualities(end)], [x_in(N2,1) x_in(N2,1)], "r--")
plot([qualities(1) qualities(end)], [x_in(N2,end) x_in(N2,end)], "r--")
plot([qualities(1) qualities(end)], [x_in(O2,1) x_in(O2,1)], "g--")
plot([qualities(1) qualities(end)], [x_in(O2,end) x_in(O2,end)], "g--")
plot([qualities(1) qualities(end)], 10*[x_in(Ar,1) x_in(Ar,1)], "b--")
plot([qualities(1) qualities(end)], 10*[x_in(Ar,end) x_in(Ar,end)], "b--")
hold off
xlabel("Vapor Mass Quality")
ylabel("Composition of Input Fraction")
legend(["Nitrogen", "Oxygen", "10*Argon"])
improvePlot

% Plot T_in, T_boil, and Q vs quality
figure(2)
plot(qualities, T_in, "-O")
hold on
plot(qualities, T_boil, "-O")
hold off
ylabel("Temperature (K)")
xlabel("Vapor Mass Quality")
legend(["Input", "Output"])
ylim([87.4 89.5])
improvePlot

figure(3)
plot(qualities, Q, "-O")
xlabel("Vapor Mass Quality")
ylabel("Feed-Specific Heat Transfer (kJ/kg)")
improvePlot

function [Q, T_boil, rl_boil, rv_boil, x_in, T_in] = Reboiler_cqP(x_boil,q_boil,P_boil)
% Returns heat transfer Q, temperature, output liquid density,
% output vapor density, input liquid composition, and temperature for a given
% output composition c, pressure P, and quality in boiler q

N = length(x_boil);
MW_O2 = 31.999; % (g/mol)
MW_N2 = 28.0134; % (g/mol)
MW_Ar = 39.948; % (g/mol)

if sum(x_boil) ~= 1
    fprintf("Input composition doe not sum to unity!! \n")
elseif q_boil < 0 || q_boil > 1
    fprintf("Input quality is outside of 0-1 bounds! \n")
end

% Assume liquid in reboiler is saturated
[T_boil, rl_boil, rv_boil, y_boil] = Fast_Bubble_cP(x_boil,P_boil);

% Find mass flow using quality of output
if N == 2 % Accomodate binary and ternary air
    MW = [MW_N2, MW_O2];
else 
    MW = [MW_N2, MW_O2, MW_Ar];
end
mdot_liq = 1-q_boil;
mdot_vap = q_boil;
ndot_liq = mdot_liq / sum(MW.*x_boil); % mass flow of liquid outut
ndot_vap = mdot_vap / sum(MW.*y_boil); % kg/(kg/kmol/s) = kmol/s  mol flow of vapor output
mdot_in = mdot_vap + mdot_liq;

% Determine input composition
x_in = ndot_liq*x_boil + ndot_vap*y_boil;
x_in = x_in/sum(x_in); % normalize values

% nput saturated Liquid
[T_in, rl_in, ~, ~] = Fast_Bubble_cP(x_in,P_boil);

% Specific enthalpy for each stream
h_vap = h_crT(y_boil, rv_boil, T_boil);
h_out = h_crT(x_boil, rl_boil, T_boil);
h_in = h_crT(x_in, rl_in, T_in);
Q = (h_out*mdot_liq + h_vap*mdot_vap - h_in*mdot_in)/1e3; % joules per kmol/s liquid out

end