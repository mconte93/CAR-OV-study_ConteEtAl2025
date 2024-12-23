%% Forward code for CAR OV module

% Clear workspace and initialize global variable
clear all
% close all  % Uncomment to close all figures

global b  % Declare global variable for b (burst size)

% Load experimental data
T_point = load('/Data/Data_extract_041122/Time_50_002.txt');  % Time points (in hours)
PBT030_Y1 = load('Data/Data_extract_041122/PBT030_50_002_Y1.txt');  % Data from experiment Y1
PBT030_Y2 = load('Data/Data_extract_041122/PBT030_50_002_Y2.txt');  % Data from experiment Y2
PBT030_Y3 = load('Data/Data_extract_041122/PBT030_50_002_Y3.txt');  % Data from experiment Y3

% Average across Y1, Y2, and Y3
PBT030_ave = mean([PBT030_Y1 PBT030_Y2 PBT030_Y3], 2);

% Plot the experimental data and error bars
subplot(2,2,3)
plot(T_point, PBT030_ave, 'k-', 'LineWidth', 0.5)  % Plot the average data in black
hold on

% Calculate error bars based on the range across Y1, Y2, and Y3
for i = 1:length(T_point)
    err(i) = max([PBT030_Y1(i), PBT030_Y2(i), PBT030_Y3(i)]) - ...
             min([PBT030_Y1(i), PBT030_Y2(i), PBT030_Y3(i)]);  % Range between max and min
end

% Add error bars to the plot
errorbar(T_point, PBT030_ave, err, 'k.', 'LineWidth', 0.1)

% Initial conditions based on the experimental data
time = T_point;  % Time vector
T_0 = PBT030_ave(1);  % Initial tumor value (from first data point)
C_0_cellNum = ((T_0 - 0.1602) / 0.0001946) / 50;  % Calculate initial CAR T cell number
C_0 = C_0_cellNum * 0.0001946 + 0.1602;  % Convert CAR T cell number back to concentration
I_0 = 0;  % Initial immune cell count
V_0 = 0.002;  % Initial virus concentration

% Set initial conditions vector
Init = [T_0; I_0; V_0; C_0];

% Load optimized parameters for the model
load('K_opt_Tonly.mat', 'K_opt')  % Load optimal K
K = K_opt;  % Tumor carrying capacity 

% Load optimized parameters for immune-virus interactions
load('alfak1k3delta_opt_T50.mat', 'delta_opt', 'k3_opt', 'k1_opt')
delta = delta_opt;  % CAR T death rate
theta4 = k3_opt;  % Expansion/Exhaustion CAR T rate

% Load parameters for immune dynamics
load('alfabetagamma_opt_OV002_c25.mat', 'gamma_opt', 'beta_opt')
gamma = gamma_opt;  % Tumor death rate
omega = 0.05;  % Virus decay rate (constant)
b = 25;  %burst size (100 PFU per cell)

% Load parameters for tumor growth and immune response
load('MeanParameterValues_theta123alfabeta.mat', 'beta_50_002_ave', 'theta1_50_002_ave', 'alfa_50_002_ave', ...
    'theta2_50_002_ave', 'theta3_50_002_ave')
beta = beta_50_002_ave;  % Tumor infection rate
alfa = alfa_50_002_ave;  % Tumor growth rate
theta1 = theta1_50_002_ave;  % CAR T kiiling rate of tumor cells
theta2 = theta2_50_002_ave;  % CAR T kiiling rate of infected tumor cells
theta3 = theta3_50_002_ave;  % CAR T infection rate 

% Pack parameters into a single vector for ODE function
par = [alfa, K, beta, theta1, gamma, theta2, b, theta3, omega, delta, theta4];

% Solve the system of ODEs
[t, y] = ode45(@(t, y) TumorCarOv(t, y, par), time, Init);

% Plot the model predictions
hold on
plot(t, y(:, 1) + y(:, 2), 'r', 'LineWidth', 2)  % Total tumor and immune cell population
plot(t, y(:, 1), 'r.-', 'LineWidth', 1)  % Tumor population
plot(t, y(:, 2), 'r--', 'LineWidth', 1)  % Immune cell population
plot(t, y(:, 3) / b, 'g', 'LineWidth', 2)  % Virus concentration (normalized by b)
plot(t, y(:, 4), 'b', 'LineWidth', 2)  % Immune cell concentration

% Set axis limits and labels
xlim([2.5 72])  % Time range
ylim([0 6])  % Population range
legend('','T(t)+I(t)', 'T(t)', 'I(t)', 'V(t)', 'C(t)', 'Interpreter', 'latex')
xlabel('Time [h]', 'Interpreter', 'latex')
ylabel('CI', 'Interpreter', 'latex')
title('E:T 1:50 - MOI: 0.002', 'Interpreter', 'latex')
axis square  % Equal scaling of axes

% Optional: Save the results for future use (uncomment to save)
% save('Forward_CAROV_50_002.mat', 't', 'y', 'theta1', 'alfa', 'beta', 'gamma', 'b', 'K', 'omega', 'theta2', 'theta3', 'theta4', 'delta')
% save('Areas_CAR50_OV002.mat', 'A_exp', 'A_exp_Y1', 'A_exp_Y2', 'A_exp_Y3', 'A_num')


% Define the system of ODEs for tumor-immune interaction model
function [dydt] = TumorCarOv(t, y, par)

% Global variable for burst size 
global MOI_cell

% Unpack parameters
alfa = par(1);   
K = par(2);      
beta = par(3);  
theta1 = par(4);
gamma = par(5);  
theta2 = par(6); 
MOI_cell = par(7);  
theta3 = par(8); 
omega = par(9);  
delta = par(10); 
theta4 = par(11); 

% Unpack state variables (T = tumor cells, I = immune cells, V = virus, C = CAR T cells)
T = y(1); 
I = y(2);  
V = y(3);  
C = y(4);  

% Define the system of differential equations
dT = alfa * T * (1 - (T + I) / K) - beta * T * V - theta1 * C * T;  % Tumor growth
dI = beta * T * V - gamma * I - theta2 * I * C;  % Immune cell dynamics
dV = gamma * MOI_cell * I - omega * V;  % Virus dynamics
dC = -delta * C + theta4 * C * (T + I) - theta3 * C * V;  % CAR T cell dynamics

% Return the derivatives
dydt = [dT; dI; dV; dC];
end
