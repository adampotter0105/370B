P = 3e5; % Operating pressure (assume no drop)
T = 1000 + 273.15; % Operating temperature

% TEST
FC_V_I(100, T, P)
FC_V_I(150, T, P)

% Current Sweep
I_range = linspace(0,205, 50);
n = length(I_range);

% Data Arrays
V_data = NaN(1,n);
Q_data = NaN(1,n);
xH2_data = NaN(1,n);
xO2_data = NaN(1,n);

ohmic_loss_data = NaN(1,n);
gdl_loss_data = NaN(1,n);
anode_loss_data = NaN(1,n);
cathode_loss_data = NaN(1,n);

% Collec Data Sweeping Across Current
for i = 1:n
    [V_data(i), Q_data(i), xH2_data(i), xO2_data(i), losses] = FC_V_I(I_range(i), T, P);
    ohmic_loss_data(i) = losses.ohmic_loss;
    gdl_loss_data(i) = losses.gdl_loss;
    anode_loss_data(i) = losses.anode_loss;
    cathode_loss_data(i) = losses.cathode_loss;
end