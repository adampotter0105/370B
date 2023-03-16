function [v_guess, Q, xH2, xO2, loss_total] = FC_V_I(I, T, P)
% Returns fuel cell voltage, heat flux, output H2 mole fraction, and output O2
% mole fraction for a given total current, temperature, and pressure

F = 96485; % C/mol, Faraday's constant
% Channel Parameters
channel_w = 10e-3; % meters
channel_h = 5e-5; % meters
channel_l = 0.5; % meters

diff_elements = 20; % number of differential button elements along channel
diff_area = channel_w*channel_l/diff_elements; % Area per button cell
cross_area = channel_w*channel_h; % channel cross-sectional area

%% Gas Flow 
xH2_0 = 0.97; % Input mole fraction
xO2_0 = 0.21; % Input mole fraction
flow_anode = 1; % m/s
lam_cath = 2; % ratio of needed air needed for reaction

gas = GRI30;
iH2  = speciesIndex(gas, 'H2');
iH2O = speciesIndex(gas, 'H2O');
iO2  = speciesIndex(gas, 'O2');
iN2 = speciesIndex(gas, 'N2');
nsp       = nSpecies(gas);
set(gas, "T", T, "P", P)

x_anode   = zeros(nsp, 1);
x_cathode = zeros(nsp, 1);
x_cathode(iO2, 1) = xO2_0; 
x_cathode(iN2, 1) = 1-xO2_0;
set(gas, 'X', x_cathode); 
rho_cathode = molarDensity(gas)/1e3;

% General Equation for mass change
% cross_area*flowrate*molar_density*mole_fraction = current * (delta_mol_X/delta_mol_e- ) / F
flow_cathode = lam_cath*I*0.25/(F*rho_cathode*cross_area*xO2_0); % Current determines how much O2 is needed
% TODO: fix this flow equation, mole fraction updates may also be incorrect

% HEAT FLUX
Q = 0; % TODO: Not sure best way to calculate this

%Iterate over elements in channel
loss_total.ohmic_loss = 0;   loss_total.cathode_loss = 0;
loss_total.anode_loss = 0;   loss_total.gdl_loss = 0;
        
% Uggghhhhh another newton raphson
tol = 1e-4;
dv = 1e-3;
max_it = 20;
it = 0;
v_guess = 0.7; % TODO: IMPROVE THIS
v_max = 1.1; % TODO: Improve This

while it < max_it
    it = it + 1;
    v_prev = v_guess;

    %Iterate through Elements
    xH2 = xH2_0;  xO2 = xO2_0; I_total = 0;
    for i = 1:diff_elements
        [i_flux, losses] = button_FC(v_guess, T, P, xH2, xO2);
        % Sum current along channel
        i_diff = i_flux*diff_area;
        I_total  = I_total + i_diff;
        loss_total.ohmic_loss = loss_total.ohmic_loss + losses.ohmic_loss;   
        loss_total.cathode_loss = loss_total.cathode_loss + losses.cathode_loss;
        loss_total.anode_loss = loss_total.anode_loss + losses.anode_loss;   
        loss_total.gdl_loss = loss_total.gdl_loss + losses.gdl_loss;
        % Account for reduced reactants and Increased Products
        % mole_fraction = (current*F*(delta_mol_X/delta_mol_e-))/(cross_area*molar_density*flowrate)
        x_anode(iH2, 1)   = xH2; x_anode(iH2O, 1)  = 1-xH2;
        x_cathode(iO2, 1) = xO2;  x_cathode(iN2, 1) = 1-xO2;
        set(gas, 'X', x_anode);   rho_anode = molarDensity(gas)/1e3;
        set(gas, 'X', x_cathode);  rho_cathode = molarDensity(gas)/1e3;

        xH2 = xH2 - i_diff*0.5/(F*cross_area*rho_anode*flow_anode);
        xO2 = xO2 - i_diff*0.5/(F*cross_area*rho_cathode*flow_cathode);
    end

    % Check for Convergence
    if abs(I_total-I)/(I+tol) < tol
        return
    end

    % Adjust with NR
    %  High
    xH2 = xH2_0;  xO2 = xO2_0; I_high = 0;
    for i = 1:diff_elements
        [i_flux, ~] = button_FC(v_guess+dv, T, P, xH2, xO2);
        i_diff = i_flux*diff_area;
        I_high  = I_high + i_diff;
        x_anode(iH2, 1)   = xH2; x_anode(iH2O, 1)  = 1-xH2;
        x_cathode(iO2, 1) = xO2;  x_cathode(iN2, 1) = 1-xO2;
        set(gas, 'X', x_anode);   rho_anode = molarDensity(gas)/1e3;
        set(gas, 'X', x_cathode);  rho_cathode = molarDensity(gas)/1e3;
        xH2 = xH2 - i_diff*0.5/(F*cross_area*rho_anode*flow_anode);
        xO2 = xO2 - i_diff*0.5/(F*cross_area*rho_cathode*flow_cathode);
    end
    % Low
    xH2 = xH2_0;  xO2 = xO2_0; I_low = 0;
    for i = 1:diff_elements
        [i_flux, ~] = button_FC(v_guess-dv, T, P, xH2, xO2);
        i_diff = i_flux*diff_area;
        I_low  = I_low + i_diff;
        x_anode(iH2, 1)   = xH2; x_anode(iH2O, 1)  = 1-xH2;
        x_cathode(iO2, 1) = xO2;  x_cathode(iN2, 1) = 1-xO2;
        set(gas, 'X', x_anode);   rho_anode = molarDensity(gas)/1e3;
        set(gas, 'X', x_cathode);  rho_cathode = molarDensity(gas)/1e3;
        xH2 = xH2 - i_diff*0.5/(F*cross_area*rho_anode*flow_anode);
        xO2 = xO2 - i_diff*0.5/(F*cross_area*rho_cathode*flow_cathode);
    end

    dedv = ((I_high-I) - (I_low-I))/(2*dv);
    i_err = I_total - I;
    v_guess = v_guess - i_err/dedv;

    % Bisection if guess out of bounds
    if v_guess < 0
        v_guess = v_prev/2;
    elseif v_guess > v_max
        v_guess = (v_max - v_prev)/2 + v_prev;
    end


end

fprintf("FC_I_V Error: Failed to converge channel voltage for a given current! \n")