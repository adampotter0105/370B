function [V, Q, xH2, xO2] = FC_V_I(I, T, P)
% Returns fuel cell voltage, heat flux, output H2 mole fraction, and output O2
% mole fraction for a given total current, temperature, and pressure

% Channel Parameters
channel_w = 10e-3; % meters
channel_h = 5e-5; % meters
channel_l = 0.5; % meters

diff_elements = 100; % number of differential button elements along channel
diff_area = channel_w*channel_l/diff_elements; % Area per button cell
diff_vol = diff_area*channel_h; % volume per button cell

% Gas Flow Params
xH2 = 0.97; xO2 = 0.21;
flow_anode = 1; % m/s
lam_cath = 2; % ratio of needed air for reaction

%Iterate over elements in channel
V_guess = 1; % IMPROVE THIS
I_total = 0;
loss_total.ohmic_loss = 0;   loss_total.cathode_loss = 0;
loss_total.anode_loss = 0;   loss_total.gdl_loss = 0;
        

%Iterate through Elements
for i = 1:diff_elements

    [i_flux, losses] = button_FC(V_guess, T, P, xH2, xO2);

    % Sum current along channel
    I_total  = I_total + i_flux*diff_area;
    loss_total = loss_total + losses;

    % Account for reduced produts and reactants
    

end