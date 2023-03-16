% Channel Parameters
channel_w = 10e-3; % meters
channel_h = 5e-5; % meters
channel_l = 0.5; % meters

diff_elements = 100; % number of differential button elements along channel
diff_area = channel_w*channel_l/diff_elements; % Area per button cell
diff_vol = diff_area*chennel_h; % volume per button cell

P = 3e5; % Operating pressure (assume no drop)
T = 1000 + 273.15; % Operating temperature

flow_anode = 1; % m/s
lam_cath = 2; % ratio of needed air for reaction

%Iterate over elements in channel
V_target = 1;
xH2 = 0.97; xO2 = 0.21;
I_total = 0;
for i = 1:diff_elements
    i_flux = button_FC(V_target, T, P, xH2, xO2);

    % Sum current along channel
    I_total  = I_total + i_flux*diff_area;

    % Account for reduced produts and reactants
    

end

