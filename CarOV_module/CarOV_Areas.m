% Close all figures and clear workspace
close all
clear all

% Declare global variable for burst size
global b

% Load experimental data for time points and tumor dynamics
T_point = load('Data/Data_Area/Time_50_00012.txt');  % Time points (in hours)
PBT030_Y1 = load('Data/Data_Area/PBT030_50_00012_Y1.txt');  % Experimental data for Y1
PBT030_Y2 = load('Data/Data_Area/PBT030_50_00012_Y2.txt');  % Experimental data for Y2
PBT030_Y3 = load('Data/Data_Area/PBT030_50_00012_Y3.txt');  % Experimental data for Y3

% Calculate the average of the three data sets (Y1, Y2, Y3)
PBT030_ave = mean([PBT030_Y1, PBT030_Y2, PBT030_Y3], 2);

% Plot the average data and error bars
plot(T_point, PBT030_ave, 'k-', 'LineWidth', 0.5)  % Plot the average curve in black
hold on

% Calculate the error bars (range between min and max values for each time point)
for i = 1:length(T_point)
    err(i) = max([PBT030_Y1(i), PBT030_Y2(i), PBT030_Y3(i)]) - ...
             min([PBT030_Y1(i), PBT030_Y2(i), PBT030_Y3(i)]);  % Max-min range
end

% Plot the error bars
errorbar(T_point, PBT030_ave, err, 'k.', 'LineWidth', 0.1)

% Set the time vector for the ODE solution
time = T_point;

% Initial conditions from the experimental data
T_0 = PBT030_ave(1);  % Initial tumor value (from the first data point)
C_0_cellNum = ((T_0 - 0.1602) / 0.0001946) / 50;  % Initial CAR T cells number
C_0 = C_0_cellNum * 0.0001946 + 0.1602;  % Convert back to concentration
I_0 = 0;  % Initial immune cell count
V_0 = 0.00012;  % Initial virus concentration

% Set the initial conditions for the ODE solver
Init = [T_0; I_0; V_0; C_0];

% Load optimal parameters for the tumor model
load('K_opt_Tonly.mat', 'K_opt')
K = K_opt;  % Tumor carrying capacity 

% Load optimized parameters
load('alfak1k3delta_opt_T50.mat', 'delta_opt', 'k3_opt', 'k1_opt', 'alfa_opt')
delta = delta_opt;  
theta4 = k3_opt;  

load('alfabetagamma_opt_OV03_c25.mat', 'gamma_opt', 'beta_opt', 'alfa_opt')
gamma_03 = gamma_opt;  %
load('alfabetagamma_opt_OV002_c25.mat', 'gamma_opt', 'beta_opt', 'alfa_opt')
gamma_002 = gamma_opt;  
load('alfabetagamma_opt_OV0008_c25.mat', 'gamma_opt', 'beta_opt', 'alfa_opt')
gamma_0008 = gamma_opt;  

% Fit a polynomial to the gamma values for different MOI
x_gamma = [0 0.0008, 0.002, 0.03];
y_gamma = [0 gamma_0008, gamma_002, gamma_03];
p_gamma = polyfit(x_gamma, y_gamma, 2);  % Fit a 2nd degree polynomial
gamma_00012 = polyval(p_gamma, 0.00012);  % Estimate gamma for MOI = 0.00012

% Other constants
omega = 0.05;  % Virus decay rate (constant)
b = 250;  % burst size (100 PFU/cell)

% Load parameter values for tumor dynamics and immune responses
load('MeanParameterValues_theta12alfabeta_theta3.mat', ...
    'beta_50_0008_ave', 'theta1_50_0008_ave', 'alfa_50_0008_ave', ...
    'theta2_50_0008_ave', 'theta3_50_0008_ave')

% Tumor model parameters for different MOIs
beta_0008 = beta_50_0008_ave;
alfa_0008 = alfa_50_0008_ave;
theta1_0008 = theta1_50_0008_ave;
theta2_0008 = theta2_50_0008_ave;
theta3_0008 = theta3_50_0008_ave;

% Load parameters for MOI = 0.002
load('MeanParameterValues_theta12alfabeta_theta3.mat', ...
    'beta_50_002_ave', 'theta1_50_002_ave', 'alfa_50_002_ave', ...
    'theta2_50_002_ave', 'theta3_50_002_ave')

beta_002 = beta_50_002_ave;
alfa_002 = alfa_50_002_ave;
theta1_002 = theta1_50_002_ave;
theta2_002 = theta2_50_002_ave;
theta3_002 = theta3_50_002_ave;

% Fit polynomials to parameters (alfa, beta, theta1, etc.)
x_alfa = [0 0.0008, 0.002];
y_alfa = [0 alfa_0008, alfa_002];
p_alfa = polyfit(x_alfa, y_alfa, 2);
alfa_00012 = polyval(p_alfa, 0.00012);  % Estimate alfa for MOI = 0.00012

x_beta = [0 0.0008, 0.002];
y_beta = [0 beta_0008, beta_002];
p_beta = polyfit(x_beta, y_beta, 2);
beta_00012 = polyval(p_beta, 0.00012);  % Estimate beta for MOI = 0.00012

x_theta1 = [0 0.0008, 0.002];
y_theta1 = [0 theta1_0008, theta1_002];
p_theta1 = polyfit(x_theta1, y_theta1, 2);
theta1_00012 = polyval(p_theta1, 0.00012);  % Estimate theta1 for MOI = 0.00012

x_theta2 = [0 0.0008, 0.002];
y_theta2 = [0 theta2_0008, theta2_002];
M = fit(x_theta2', y_theta2', 'poly2', 'lower', [0, 0, 0]);
theta2_00012 = feval(M, 0.00012);  % Estimate theta2 for MOI = 0.00012

x_theta3 = [0 0.0008, 0.002];
y_theta3 = [0 theta3_0008, theta3_002];
p_theta3 = polyfit(x_theta3, y_theta3, 2);
theta3_00012 = polyval(p_theta3, 0.00012);  % Estimate theta3 for MOI = 0.00012

% Set all model parameters for the ODE solver
par = [alfa_00012, K, beta_00012, theta1_00012, gamma_00012, theta2_00012, b, theta3_00012, omega, delta, theta4];

% Solve the system of ODEs using ode45
[t, y] = ode45(@(t, y) TumorCarOv(t, y, par), time, Init);

% Plot the model results
hold on
plot(t, y(:, 1) + y(:, 2), 'r', 'LineWidth', 2)  % Tumor + Immune cell population
plot(t, y(:, 1), 'r.-', 'LineWidth', 1)  % Tumor population
plot(t, y(:, 2), 'r--', 'LineWidth', 1)  % Immune cell population
plot(t, y(:, 3) / b, 'g', 'LineWidth', 2)  % Virus concentration (normalized)
plot(t, y(:, 4), 'b', 'LineWidth', 2)  %CAR T population

% Set plot limits and labels
xlim([2.5 72])  % Time range
ylim([0 6])  % Population range
legend('', '', 'T(t)+I(t)', 'T(t)', 'I(t)', 'V(t)', 'C(t)', 'Interpreter', 'latex')
xlabel('Time [h]', 'Interpreter', 'latex')
ylabel('CI', 'Interpreter', 'latex')
title('E:T 1:50 - $V_0$=0.00012 MOI', 'Interpreter', 'latex')
axis square  % Equal scaling of axes

% Define the system of ODEs for tumor-immune interaction model
function [dydt] = TumorCarOv(t, y, par)

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

    % Unpack state variables
    T = y(1);  % Tumor population
    I = y(2);  % Immune population
    V = y(3);  % Virus population
    C = y(4);  % CAR T population

    % Differential equations
    dT = alfa * T * (1 - (T + I) / K) - beta * T * V - theta1 * C * T;  % Tumor growth dynamics
    dI = beta * T * V - gamma * I - theta2 * I * C;  % Immune cell dynamics
    dV = gamma * MOI_cell * I - omega * V;  % Virus dynamics
    dC = -delta * C + theta4 * C * (T + I) - theta3 * C * V;  % CAR T dynamics

    % Return derivatives
    dydt = [dT; dI; dV; dC];
end

