%% Code for forward OV module

clear all;  % Clear workspace variables

% Load time and experimental data for tumor dynamics
T_point = load('Data/Data_extract_041122/Time_OV0008.txt');  % Time points (hours)
PBT030_Y1 = load('Data/Data_extract_041122/PBT030_OV0008_Y1.txt');  % Data from set 1
PBT030_Y2 = load('Data/Data_extract_041122/PBT030_OV0008_Y2.txt');  % Data from set 2
PBT030_Y3 = load('Data/Data_extract_041122/PBT030_OV0008_Y3.txt');  % Data from set 3

% Average experimental data from all three sets
PBT030_ave = mean([PBT030_Y1, PBT030_Y2, PBT030_Y3], 2);

% Plot the experimental data with error bars
figure(1)
subplot(2, 2, 3)
plot(T_point, PBT030_ave, 'k-', 'LineWidth', 0.5);  % Plot the average data
hold on;

% Calculate error bars (difference between max and min for each time point)
for i = 1:length(T_point)
    err(i) = max([PBT030_Y1(i), PBT030_Y2(i), PBT030_Y3(i)]) - min([PBT030_Y1(i), PBT030_Y2(i), PBT030_Y3(i)]);
end

% Plot error bars for the experimental data
errorbar(T_point, PBT030_ave, err, 'k.', 'LineWidth', 0.1);

% Load previously optimized parameters
load('K_opt_Tonly.mat');  % Load tumor carrying capacity (K)
K = K_opt;  % Assign value of K

load('Extended_Burst_size_analysis/alfabetagamma_opt_OV0008_c25.mat');  % Load optimized parameters for alfa, beta, gamma
alfa = alfa_opt;
gamma = gamma_opt;
beta = beta_opt;

b = 25;  % Burst size (MOI_ICI)
omega = 0.05;  % Virus clearance rate (1/h)

% Initial conditions for the ODE model
time = T_point;
T_0 = PBT030_ave(1);  % Initial tumor concentration (first point of the average data)
I_0 = 0;  % Initial immune cell concentration (assumed to be zero)
V_0 = 0.0008;  % Initial virus concentration (example value)
Init = [T_0; I_0; V_0];  % Initial conditions vector

% Solve the system of ODEs using ode45
[t, y] = ode45(@(t, y) Ov(t, y, [alfa, K, beta, gamma, b, omega]), time, Init, []);

% Plot results of the model
figure(1)
subplot(2, 2, 3)
hold on;
plot(t, y(:, 1) + y(:, 2), 'r', 'LineWidth', 2);  % Tumor + Non-infected cells over time
plot(t, y(:, 1), 'b.-', 'LineWidth', 0.5);  % Tumor cells over time
plot(t, y(:, 2), 'b--', 'LineWidth', 1);  % Non-infected immune cells over time
plot(t, y(:, 3) / b, 'g', 'LineWidth', 2);  % Oncolytic virus over time (scaled by burst size)

% Customize the plot
xlim([T_point(1) 72]);
ylim([0 5]);
yticks([0 1 2 3 4 5]);
xticks([12 24 36 48 60 72]);
legend('Tumor + Immune cells', 'Tumor cells', 'Non-Infected cells', 'Infected cells', 'Oncolytic Virus', ...
       'Interpreter', 'latex', 'FontSize', 11, 'Location', 'northwest');
xlabel('Time [h]', 'Interpreter', 'latex', 'FontSize', 13);
ylabel('[CI]', 'Interpreter', 'latex', 'FontSize', 13);
title('MOI = 0.0008', 'Interpreter', 'latex', 'FontSize', 13);
axis square;

% Uncomment the following lines to save data
% save('Forward_OV_0008.mat', 't', 'y', 'alfa', 'beta', 'gamma', 'b', 'K', 'omega');

%%% ODE function for Tumor-OV dynamics
function [dydt] = Ov(t, y, par)

% Unpack parameters
alfa = par(1);  % Tumor growth rate
K = par(2);     % Tumor carrying capacity
beta = par(3);  % Infection rate
gamma = par(4); % Immune cell response rate
b = par(5);     % Burst size (MOI_ICI)
omega = par(6); % Virus clearance rate

% State variables
T = y(1);  % Tumor cell concentration
I = y(2);  % Infected immune cell concentration
V = y(3);  % Virus concentration

% Define the system of differential equations
dTdt = alfa * T * (1 - (T + I) / K) - beta * T * V;  % Tumor dynamics
dIdt = beta * T * V - gamma * I;  % Immune cell dynamics
dVdt = b * gamma * I - omega * V;  % Virus dynamics

% Return the derivatives
dydt = [dTdt; dIdt; dVdt];

end