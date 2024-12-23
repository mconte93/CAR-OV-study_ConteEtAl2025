%% Forward code for CAR module

clear all;
% close all; % Uncomment to close all figures before running (if needed)

% Load experimental data (in hours)
T_point = load('/Data/Data_extract_041122/Time_carT10.txt'); % Time points
PBT030_Y1 = load('/Data/Data_extract_041122/PBT030_carT10_Y1.txt'); % PBT030 dataset for Year 1
PBT030_Y2 = load('/Data/Data_extract_041122/PBT030_carT10_Y2.txt'); % PBT030 dataset for Year 2
PBT030_Y3 = load('/Data/Data_extract_041122/PBT030_carT10_Y3.txt'); % PBT030 dataset for Year 3

% Calculate average PBT030 across the three replicates
PBT030_ave = mean([PBT030_Y1, PBT030_Y2, PBT030_Y3], 2);

% Plot experimental data with error bars
figure(1)
subplot(2,2,3)
plot(T_point, PBT030_ave, 'k-', 'LineWidth', 0.5); % Plot average PBT030
hold on;

% Calculate error for each data point as the range (max - min) of the replicates 
err = zeros(size(T_point));
for i = 1:length(T_point)
    err(i) = max([PBT030_Y1(i), PBT030_Y2(i), PBT030_Y3(i)]) - min([PBT030_Y1(i), PBT030_Y2(i), PBT030_Y3(i)]);
end

% Plot error bars
errorbar(T_point, PBT030_ave, err, 'k.', 'LineWidth', 0.1);

% Load model parameters from .mat files
load('MeanParameterValues_alfadeltak1k3.mat');  % Contains model parameters like alfa, delta, theta_TC, mu
load('K_opt_Tonly.mat'); % Contains optimal K value

% Assign parameter values to variables for clarity
alfa = alfa_10_ave;
K = K_opt;
k1 = k1_10_ave;  %\theta_TC
k3 = k3_10_ave;  %\mu
delta = delta_10_ave;

% Set initial conditions and time span for ODE solver
T_0 = PBT030_ave(1);  % Initial tumor cell count
C_0_cellNum = ((PBT030_ave(1) - 0.1602) / 0.0001946) / 10; % Cell number from concentration
C_0 = C_0_cellNum * 0.0001946 + 0.1602;  % Initial Car T-cell concentration

Init = [T_0; C_0];  % Initial conditions vector: [Tumor cells, Car T-cells]

% Solve ODEs using ode45
[t, y] = ode45(@(t, y) CarOv(t, y, [alfa, K, k1, k3, delta]), T_point, Init);

% Plot the results (Tumor and Car T-cells over time)
subplot(2,2,3);
plot(t, y(:, 1), 'r', 'LineWidth', 2); % Tumor cells (T) over time
plot(t, y(:, 2), 'b', 'LineWidth', 2); % Car T-cells (C) over time

% Customize plot appearance
xlim([T_point(1) 72]);
ylim([0 5]);
yticks([0 1 2 3 4 5]);
xticks([12 24 36 48 60 72]);

% Add legend, labels, and title
legend('Tumor cells', 'Car T-cells', 'Interpreter', 'latex', 'FontSize', 13);
xlabel('Time [h]', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('[ICI]', 'Interpreter', 'latex', 'FontSize', 13);
title('E:T ratio - 1:10', 'Interpreter', 'latex', 'FontSize', 13);

% Set axis square for better aspect ratio
axis square;

% Uncomment to save results if needed
% save('Forward_CAR_10.mat', 't', 'y', 'alfa', 'k1', 'k3', 'delta', 'K');

% Define the system of ODEs for Tumor and Car T-cell dynamics
function dydt = CarOv(t, y, par)
    alfa = par(1);
    K = par(2);
    k1 = par(3);
    k3 = par(4);
    delta = par(5);

    T = y(1);  % Tumor cells
    C = y(2);  % Car T-cells

    % Differential equations
    dTdt = alfa * T * (1 - T / K) - k1 * T * C;  % Tumor cell growth and interaction with Car T-cells
    dCdt = -delta * C + k3 * T * C;  % Car T-cell decay and expansion/exhaustion by tumor cells

    dydt = [dTdt; dCdt];  % Return the derivatives
end

