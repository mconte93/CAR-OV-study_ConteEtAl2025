%%% Estimation of beta parameter only

close all;
clear all;

% Define global variables for the model parameters
global K MOI_ICI omega alfa gamma

% Load experimental data for tumor-only condition from 2.5h to 72h
T_point = load('Data/Data_extract_041122/Time_OV03.txt');  % Time points (in hours)
PBT030_Y1 = load('Data/Data_extract_041122/PBT030_OV03_Y1.txt');  % Tumor data Y1
PBT030_Y2 = load('Data/Data_extract_041122/PBT030_OV03_Y2.txt');  % Tumor data Y2
PBT030_Y3 = load('Data/Data_extract_041122/PBT030_OV03_Y3.txt');  % Tumor data Y3

% Calculate the average tumor data across three measurements (Y1, Y2, Y3)
PBT030_ave = mean([PBT030_Y1 PBT030_Y2 PBT030_Y3], 2);

% Plot the experimental tumor data with error bars
figure(2);
plot(T_point, PBT030_ave, 'k-', 'LineWidth', 0.5);
hold on;
for i = 1:size(T_point)
    err(i) = max([PBT030_Y1(i), PBT030_Y2(i), PBT030_Y3(i)]) - min([PBT030_Y1(i), PBT030_Y2(i), PBT030_Y3(i)]);
end
errorbar(T_point, PBT030_ave, err, 'k.', 'LineWidth', 0.1);

% Set initial conditions for the ODE solver
time = T_point;
T_0 = PBT030_ave(1);
I_0 = 0;
V_0 = 0.03;
Init = [T_0; I_0; V_0];

% Load the optimized tumor carrying capacity (K)
load('K_opt_Tonly.mat', 'K_opt');
K = K_opt;  % Tumor carrying capacity

% Set model parameters
MOI_ICI = 25;  % Burst size
omega = 0.05;    % Virus clearance rate (1/h) from Hillen et al. 2023

% Load the average parameter values for alfa and gamma
load('Mean_alfabetagamma_002.mat', 'alfa_ave_002', 'gamma_ave_002');
load('Mean_alfabetagamma_0008.mat', 'alfa_ave_0008', 'gamma_ave_0008');
load('Mean_alfabetagamma_03.mat', 'alfa_ave_03', 'gamma_ave_03');

% Average the parameters across different experimental datasets
alfa = mean([alfa_ave_002, alfa_ave_0008, alfa_ave_03]);
gamma = mean([gamma_ave_002, gamma_ave_0008, gamma_ave_03]);

% Initial guess for the beta parameter (infection rate)
beta_0 = 0.1;
x0 = beta_0;

% Set options for the particle swarm optimization
options = optimoptions('particleswarm', 'SwarmSize', 200, 'FunctionTolerance', 1e-10, ...
    'HybridFcn', @fmincon, 'MaxIterations', 500, 'ObjectiveLimit', 1e-10, ...
    'InitialSwarmMatrix', x0, 'PlotFcn', @pswplotbestf);

% Set lower and upper bounds for beta parameter
lb = 0;
ub = 10;
nvar = 1;  % Number of variables (just beta)

% Define the objective function for parameter estimation
fun = @fun_beta_PS;

% Run the particle swarm optimization to estimate the optimal beta parameter
[x, fval, exitflag, output] = particleswarm(fun, nvar, lb, ub, options);

% Display the results of the optimization
fprintf('particleswarm reached the value %f using %d function evaluations.\n', fval, output.funccount);

% Extract the optimal beta parameter from the optimization result
beta_opt = x(1);

% Solve the ODEs using the optimal beta parameter
[t, y] = ode23(@(t, y) TumorOV(t, y, [alfa, K, beta_opt, gamma, omega]), time, Init, []);
hold on;

% Plot the model simulation results
figure(2);
plot(t, y(:, 1) + y(:, 2), 'r', t, y(:, 3) / MOI_ICI, 'b', 'LineWidth', 2);

% Set plot limits and labels
xlim([0, 72]);
ylim([0, 6]);
legend('Tumor cells + Infected cells', 'Oncolytic Virus (normalized)', 'Location', 'best');
xlabel('Time (hours)');
ylabel('ICIs');
title('MOI 0.03');
axis square;

% Save the optimal beta parameter
save('beta_opt_OV03_c25.mat', 'beta_opt');

% Objective function for particle swarm optimization
function SSE = fun_beta_PS(x)

global K omega gamma alfa

% Load experimental data
T_point = load('Data/Data_extract_041122/Time_OV03.txt');
PBT030_Y1 = load('Data/Data_extract_041122/PBT030_OV03_Y1.txt');
PBT030_Y2 = load('Data/Data_extract_041122/PBT030_OV03_Y2.txt');
PBT030_Y3 = load('Data/Data_extract_041122/PBT030_OV03_Y3.txt');
PBT030_ave = mean([PBT030_Y1, PBT030_Y2, PBT030_Y3], 2);

xdata = T_point;
ydata = PBT030_ave;

% Extract beta from the optimization variable
beta = x(1);

% Set parameters for the model
par = [alfa, K, beta, gamma, omega];
T_0 = PBT030_ave(1);
I_0 = 0;
V_0 = 0.03;
Init = [T_0; I_0; V_0];

% Solve the model ODEs
[t, y] = ode45(@(t, y) TumorOV(t, y, par), xdata, Init, []);

% Calculate the sum of squared errors (SSE)
SSE = sum((ydata - (y(:, 1) + y(:, 2))).^2);

end

% Tumor-immune-virus model ODEs
function [dydt] = TumorOV(t, y, par)

global MOI_ICI

% Extract model parameters
alfa = par(1);
K = par(2);
beta = par(3);
gamma = par(4);
omega = par(5);

% Extract state variables: tumor cells (T), infected cells (I), virus (V)
T = y(1);
I = y(2);
V = y(3);

% Tumor growth and infection dynamics
dTdt = alfa * T * (1 - (T + I) / K) - beta * T * V;  % Tumor cell dynamics
dIdt = beta * T * V - gamma * I;  % Infected cell dynamics
dVdt = MOI_ICI * gamma * I - omega * V;  % Virus dynamics

% Return the derivatives
dydt = [dTdt; dIdt; dVdt];

end