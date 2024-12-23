%% Parameter Estimation - Tumor Growth and CAR T - OV Model

close all
clear all

% Define global variables
global K theta4 delta gamma MOI_cell omega

%%% Load Tumor Data (from 2.5h to 72h)
T_point = load('Data/Data_extract_041122/Time_50_0008.txt');  % Time points (in hours)
PBT030_Y1 = load('Data/Data_extract_041122/PBT030_50_0008_Y1.txt');  % Experimental data from Group 1 (hours)
PBT030_Y2 = load('Data/Data_extract_041122/PBT030_50_0008_Y2.txt');  % Experimental data from Group 2 (hours)
PBT030_Y3 = load('Data/Data_extract_041122/PBT030_50_0008_Y3.txt');  % Experimental data from Group 3 (hours)

% Compute the average of the three experimental groups
PBT030_ave = mean([PBT030_Y1, PBT030_Y2, PBT030_Y3], 2);

% Plot the average tumor growth data with error bars
figure(2)
plot(T_point, PBT030_ave, 'k-', 'LineWidth', 0.5)
hold on

% Calculate error bars (difference between max and min across the 3 groups)
for i = 1:size(T_point)
    err(i) = max([PBT030_Y1(i), PBT030_Y2(i), PBT030_Y3(i)]) - min([PBT030_Y1(i), PBT030_Y2(i), PBT030_Y3(i)]);
end
errorbar(T_point, PBT030_ave, err, 'k.', 'LineWidth', 0.1)

% Initial conditions for the ODE system
time = T_point;  % Time points for simulation
T_0 = PBT030_ave(1);  % Initial tumor size 
C_0_cellNum = ((T_0 - 0.1602) / 0.0001946) / 50;  % Calculate initial CAR T cell number based on a specific conversion factor
C_0 = C_0_cellNum * 0.0001946 + 0.1602;  % Adjust initial CAR T cell concentration
I_0 = 0;  % Initial infected cells
V_0 = 0.0008;  % Initial viral load

% Initial state vector for the ODE solver: [Tumor, Infected, Virus, CAR-T Cells]
Init = [T_0; I_0; V_0; C_0];

% Load pre-optimized parameters
load('K_opt_Tonly.mat', 'K_opt')
K = K_opt;  % Carrying capacity for tumor 

load('alfak1k3delta_opt_T50.mat', 'delta_opt', 'k3_opt', 'k1_opt')
delta = delta_opt;  % Tumor death rate due to CAR T-cells interaction
theta4 = k3_opt;  % Interaction rate parameter (mu)

% Load virus infection parameters
load('alfabetagamma_opt_OV0008_c25.mat', 'gamma_opt', 'beta_opt')
gamma = gamma_opt;  % tumor death rate
omega = 0.05;  % Virus clearance rate
MOI_cell = 25;  % Burst size

% Set initial guesses for the parameters to be optimized
theta1_0 = 0.01;  % CAR-T cells killing rate of tumor cells
theta2_0 = 0.01;  % CAR-T cells killing rate of infected cells
theta3_0 = 0.01;  % Infection rate of CAR-T cells
alfa_0 = 0.1;     % Tumor growth rate
beta_0 = 0.01;    % Virus infection rate of tumor cells

x0 = [theta1_0, theta2_0, theta3_0, alfa_0, beta_0];  % Initial guess vector

% Set up Particle Swarm Optimization (PSO) options
options = optimoptions('particleswarm', ...
    'SwarmSize', 250, 'FunctionTolerance', 1e-6, 'HybridFcn', @fmincon, ...
    'MaxIterations', 200, 'ObjectiveLimit', 1e-6, 'InitialSwarmMatrix', x0, ...
    'PlotFcn', @pswplotbestf);

% Define lower and upper bounds for the parameters
lb = [0 0 0 0 0];  % Lower bounds for parameters
ub = [1 0.999 1 1 1];  % Upper bounds for parameters

nvar = 5;  % Number of parameters to optimize

% Define the objective function for optimization
fun = @fun_theta123alfabeta_PS;

% Run the Particle Swarm Optimization
[x, fval, exitflag, output] = particleswarm(fun, nvar, lb, ub, options);

% Display the optimization results
formatstring = 'particleswarm reached the value %f using %d function evaluations.\n';
fprintf(formatstring, fval, output.funccount)

% Extract the optimal parameter values
theta1_opt = x(1);
theta2_opt = x(2);
theta3_opt = x(3);
alfa_opt = x(4);
beta_opt = x(5);

% Store optimized parameters for future use
par = [alfa_opt, K, beta_opt, theta1_opt, gamma, theta2_opt, MOI_cell, theta3_opt, omega, delta, theta4];

% Solve the ODE system with optimized parameters
[t, y] = ode45(@(t, y) TumorCarOv(t, y, par), time, Init, []);

% Plot the simulation results
figure(2)
hold on
plot(t, y(:,1) + y(:,2), 'r', 'LineWidth', 2)  % Tumor + infected cells
plot(t, y(:,3) / MOI_cell, 'b', 'LineWidth', 2)  % Virus normalized by MOI
plot(t, y(:,4), 'g', 'LineWidth', 2)  % CAR T-cells

% Customize the plot
xlim([0 72])
ylim([0 6])
legend('Tumor + Infected', 'Virus', 'CAR T-cells')
xlabel('Time (hours)')
ylabel('ICI/MOI')
title('E:T Ratio - 1:50, MOI: 0.0008')

% Save the optimized parameters to a file
save('theta123alfabeta_50_0008.mat', 'theta1_opt', 'theta2_opt', 'theta3_opt', 'alfa_opt', 'beta_opt')

%% Objective function for optimization - Sum of Squared Errors (SSE)

function SSE = fun_theta123alfabeta_PS(x)

    global K theta4 delta gamma MOI_cell omega

    % Load experimental data again
    T_point = load('Data/Data_extract_041122/Time_50_0008.txt');  % Time points
    PBT030_Y3 = load('Data/Data_extract_041122/PBT030_50_0008_Y3.txt');  % Group 3 experimental data

    PBT030_ave = mean([PBT030_Y1, PBT030_Y2, PBT030_Y3], 2);  % Mean of experimental data
    xdata = T_point;
    ydata = PBT030_ave;  % Experimental data to fit the model

    % Initial tumor size (from experimental data)
    T_0 = PBT030_ave(1);

    % Extract parameters from optimization vector
    theta1 = x(1);
    theta2 = x(2);
    theta3 = x(3);
    alfa = x(4);
    beta = x(5);

    % Set the model parameters
    par = [alfa, K, beta, theta1, gamma, theta2, MOI_cell, theta3, omega, delta, theta4];

    % Calculate initial conditions for ODE
    C_0_cellNum = ((T_0 - 0.1602) / 0.0001946) / 50;
    C_0 = C_0_cellNum * 0.0001946 + 0.1602;
    I_0 = 0;
    V_0 = 0.0008;

    % Initial state vector for the ODE solver
    Init = [T_0; I_0; V_0; C_0];

    % Solve the ODE system
    [t, y] = ode45(@(t, y) TumorCarOv(t, y, par), xdata, Init, []);

    % Compute the sum of squared errors (SSE)
    R1 = sum((ydata - (y(:,1) + y(:,2))).^2);  % Tumor and infected cells

    % Calculate cell number at 60 hours (for comparison)
    CellNumT_60h = 4 * (ydata(220) - 0.1602) / 0.0001946;
    CellInC_60h = CellNumT_60h * 0.0001946 + 0.1602;
    R2 = (CellInC_60h - y(223, 4))^2;  % CAR T-cell comparison at 60 hours

    % Total SSE
    SSE = R1 + R2;
end

%% Tumor-Car T-cell-OV ODE System

function [dydt] = TumorCarOv(t, y, par)

    global K theta4 delta gamma MOI_cell omega

    % Extract parameters
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

    % Extract state variables
    T = y(1);  % Tumor cells
    I = y(2);  % Infected tumor cells
    V = y(3);  % Virus particles
    C = y(4);  % CAR T-cells

    % ODE system
    dT = alfa * T * (1 - (T + I) / K) - beta * T * V - theta1 * C * T;  % Tumor growth rate
    dI = beta * T * V - gamma * I - theta2 * I * C;  % Infected tumor cell dynamics
    dV = gamma * MOI_cell * I - omega * V;  % Virus dynamics
    dC = -delta * C - theta3 * C * V + theta4 * C * (T + I);  % CAR T-cell dynamics

    % Return derivatives
    dydt = [dT; dI; dV; dC];
end
