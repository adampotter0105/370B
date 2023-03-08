function [Q, vapout, liqin] = Reboiler_xqP(x_boil,q_boil,P_boil)
% Returns heat transfer Q, temperature, output liquid density,
% output vapor density, input liquid composition, and temperature for a given
% output composition c, pressure P, and quality in boiler q

% Inputs: vapin: P, mdot, ndot, h, c (1xN array)
%           liqout: P, mdot, ndot, h, c (1xN array)
% Outputs: vapout: P, mdot, ndot, h, c (1xN array)
%           liqin: P, mdot, ndot, h, c (1xN array)

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
ndot_in = ndot_vap + ndot_liq;

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

% Package into structs
vapout.P = P_boil;
vapout.h = h_vap;
vapout.c = y_boil;
vapout.mdot = mdot_vap;
vapout.ndot = ndot_vap;

liqin.P = P_boil;
liqin.h = h_in;
liqin.c = x_in;
liqin.mdot = mdot_in;
liqin.ndot = ndot_in;


end