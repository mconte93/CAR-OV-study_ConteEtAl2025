%% Parameter Estimation for OV Module

clear all;
close all;

% Global variables for parameters
global K MOI_ICI omega;

% Load tumor data (PBT030) from experimental results for different time points
T_point = load('Data/Data_extract_041122/Time_OV0008.txt');  % Time points (in hours)
PBT030_Y1 = load('Data/Data_extract_041122/PBT030_OV0008_Y1.txt');  % Data from experimental set 1
PBT030_Y2 = load('Data/Data_extract_041122/PBT030_OV0008_Y2.txt');  % Data from experimental set 2
PBT030_Y3 = load('Data/Data_extract_041122/PBT030_OV0008_Y3.txt');  % Data from experimental set 3

% Compute the average of the three data sets
PBT030_ave = mean([PBT030_Y1, PBT030_Y2, PBT030_Y3], 2);

% Plot the average tumor data
figure(2)
plot(T_point, PBT030_ave, 'k-', 'LineWidth', 0.5);
hold on;

% Calculate the error bars (difference between max and min of each data set)
for i = 1:length(T_point)
    err(i) = max([PBT030_Y1(i), PBT030_Y2(i), PBT030_Y3(i)]) - min([PBT030_Y1(i), PBT030_Y2(i), PBT030_Y3(i)]);
end

% Plot the data with error bars
errorbar(T_point, PBT030_ave, err, 'k.', 'LineWidth', 0.1);

% Initial conditions: Tumor, Immune cells, and Virus
time = T_point;
T_0 = PBT030_ave(1);  % Initial tumor concentration (first point of average data)
I_0 = 0;  % Initial immune cells (assumed to be zero initially)
V_0 = 0.0008;  % Initial virus concentration (as an example)

% Initial values for parameter optimization
Init = [T_0; I_0; V_0];  % Initial conditions vector

% Load the optimal carrying capacity for tumor (K) from a previous computation
load('K_opt_Tonly.mat', 'K_opt');
K = K_opt;  % Tumor carrying capacity

% Define constants for the model
MOI_ICI = 25;  % Burst size of the virus (Virus production per infected cell)
omega = 0.05;  % Virus clearance rate (1/h) from Hillen et al., 2023

% Initial guesses for parameters alfa, beta, gamma
alfa_0 = 0.15;
beta_0 = 0.1;  
gamma_0 = 0.1;

% Combine initial guesses into a single vector for optimization
x0 = [alfa_0, beta_0, gamma_0];

% Set up optimization options for particle swarm optimization (PSO)
options = optimoptions('particleswarm', 'SwarmSize', 200, 'FunctionTolerance', 1e-10, ...
                       'HybridFcn', @fmincon, 'MaxIterations', 500, 'ObjectiveLimit', 1e-10, ...
                       'InitialSwarmMatrix', x0, 'PlotFcn', @pswplotbestf);

% Set bounds for the parameters
lb = [0, 0, 0];  % Lower bounds for alfa, beta, gamma
ub = [1, 5.5, 0.5];  % Upper bounds for alfa, beta, gamma
nvar = 3;  % Number of variables (alfa, beta, gamma)

% Define the objective function for optimization (SSE calculation)
fun = @fun_alfabetagamma_PS;

% Run particle swarm optimization to find the optimal parameters
[x, fval, exitflag, output] = particleswarm(fun, nvar, lb, ub, options);

% Output results from PSO optimization
fprintf('particleswarm reached the value %f using %d function evaluations.\n', fval, output.funccount);

% Extract optimized parameters
alfa_opt = x(1);
beta_opt = x(2);
gamma_opt = x(3);

% Solve the ODE system with the optimized parameters
[t, y] = ode23(@(t, y) TumorOV(t, y, [alfa_opt, K, beta_opt, gamma_opt, omega]), time, Init);

% Plot the results
hold on;
figure(2);
plot(t, y(:, 1) + y(:, 2), 'r', 'LineWidth', 2);  % Tumor + Immune cells over time
plot(t, y(:, 3) ./ MOI_ICI, 'b', 'LineWidth', 2);  % Virus concentration over time

% Customize the plot
xlim([0, 72]);
ylim([0, 6]);
legend('Tumor cells + Immune cells', 'Virus (scaled by MOI_ICI)');
xlabel('Time (hours)');
ylabel('Concentration');
title('MOI = 0.03');
axis square;

% Save the optimized parameters
save('alfabetagamma_opt_OV0008_c25.mat', 'alfa_opt', 'beta_opt', 'gamma_opt');

%%% Objective function for PSO (calculating sum of squared errors)
function SSE = fun_alfabetagamma_PS(x)

global K omega;

% Load the data again
T_point = load('Data/Data_extract_041122/Time_OV0008.txt');
PBT030_Y1 = load('Data/Data_extract_041122/PBT030_OV0008_Y1.txt');
PBT030_Y2 = load('Data/Data_extract_041122/PBT030_OV0008_Y2.txt');
PBT030_Y3 = load('Data/Data_extract_041122/PBT030_OV0008_Y3.txt');
PBT030_ave = mean([PBT030_Y1, PBT030_Y2, PBT030_Y3], 2);

% Experimental data for fitting
xdata = T_point;
ydata = PBT030_ave;

% Extract parameters for the model
alfa = x(1); 
beta = x(2);
gamma = x(3);

% Pack parameters for passing to the ODE function
par = [alfa, K, beta, gamma, omega];

% Initial conditions (tumor, immune cells, virus)
T_0 = PBT030_ave(1);  % Initial tumor concentration
I_0 = 0;  % Initial immune cells (zero)
V_0 = 0.0008;  % Initial virus concentration

Init = [T_0; I_0; V_0];  % Initial conditions vector

% Solve the ODE system with the current parameters
[t, y] = ode45(@(t, y) TumorOV(t, y, par), xdata, Init);

% Calculate the sum of squared errors (SSE) between the model and experimental data
SSE = sum((ydata - (y(:, 1) + y(:, 2))).^2);

end

%%% ODE system for Tumor-OV dynamics
function [dydt] = TumorOV(t, y, par)

global MOI_ICI;

% Extract parameters
alfa = par(1);
K = par(2);
beta = par(3);
gamma = par(4);
omega = par(5);

% State variables: Tumor (T), Immune cells (I), Virus (V)
T = y(1);
I = y(2);
V = y(3);

% ODEs for the system
dTdt = alfa * T * (1 - (T + I) / K) - beta * T * V;  % Tumor dynamics
dIdt = beta * T * V - gamma * I;  % Immune cell dynamics
dVdt = MOI_ICI * gamma * I - omega * V;  % Virus dynamics

% Return the derivatives
dydt = [dTdt; dIdt; dVdt];

end