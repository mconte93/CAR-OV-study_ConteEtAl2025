%% Parameter estimation from fitting for Tumor and Car T-cells dynamics

close all;
clear all;

global K;

% Load experimental data (Tumor only from 2.5h to 72h)
T_point = load('Data/Data_extract_041122/Time_carT10.txt'); % Time points in hours
PBT030_Y1 = load('Data/Data_extract_041122/PBT030_carT10_Y1.txt'); % Data from replicate 1
PBT030_Y2 = load('Data/Data_extract_041122/PBT030_carT10_Y2.txt'); % Data from replicate 2
PBT030_Y3 = load('Data/Data_extract_041122/PBT030_carT10_Y3.txt'); % Data from replicate 3

% Average PBT030 data over all replicates
PBT030_ave = mean([PBT030_Y1, PBT030_Y2, PBT030_Y3], 2);

% Plot the average experimental data with error bars
figure
plot(T_point, PBT030_ave, 'k-', 'LineWidth', 0.5);
hold on;

% Calculate the error (max - min) across the three replicates for each data point
err = max([PBT030_Y1, PBT030_Y2, PBT030_Y3], [], 2) - min([PBT030_Y1, PBT030_Y2, PBT030_Y3], [], 2);
errorbar(T_point, PBT030_ave, err, 'k.', 'LineWidth', 0.1);

% Initial conditions: Tumor and Car T-cell concentrations at time 0
T_0 = PBT030_Y3(1);
C_0_cellNum = ((PBT030_Y3(1) - 0.1602) / 0.0001946) / 10;  % Convert to cell number
C_0 = C_0_cellNum * 0.0001946 + 0.1602;  % Car T-cell concentration
Init = [T_0; C_0];  % Initial conditions vector: [Tumor cells, Car T-cells]

% Load optimal K value
load('K_opt_Tonly.mat', 'K_opt');
K = K_opt;  % Tumor carrying capacity 

% Initial guesses for parameter optimization (alfa, theta_TC, mu, delta)
x0 = [0.15, 0.1, 0.1, 0.1];  % Initial guesses for [alfa, theta_TC, mu, delta]

% Optimization options for Particle Swarm
options = optimoptions('particleswarm', 'SwarmSize', 10, 'FunctionTolerance', 1e-6, 'HybridFcn', @fmincon,'MaxIterations', 500, 'ObjectiveLimit', 1e-6, 'InitialSwarmMatrix', x0);

% Lower and upper bounds for the parameters
lb = [0, 0, -5, 0];  % Lower bounds for [alfa, theta_TC, mu, delta]
ub = [5, 10, 5, 1];  % Upper bounds for [alfa, theta_TC, mu, delta]
nvar = 4;  % Number of variables to optimize

% Define the objective function for Particle Swarm Optimization
fun = @fun_alfak1k3delta_PS;

% Run the Particle Swarm Optimization to estimate parameters
[x, fval, exitflag, output] = particleswarm(fun, nvar, lb, ub, options);

% Extract optimized parameters
alfa_opt = x(1);
k1_opt = x(2);
k3_opt = x(3);
delta_opt = x(4);

% Solve the ODEs with the optimized parameters
[t, y] = ode45(@(t, y) CarOv(t, y, [x(1), K, x(2), x(3), x(4)]), T_point, Init);

% Plot the results (Tumor and Car T-cells over time)
hold on;
plot(t, y(:, 1), 'r', 'LineWidth', 2); % Tumor cells (T) over time
plot(t, y(:, 2), 'b', 'LineWidth', 2); % Car T-cells (C) over time

% Customize plot appearance
xlim([0 72]);
ylim([0 6]);
legend('Tumor cells', 'Car T-cells');
xlabel('Time [h]');
ylabel('ICI');
title('E:T ratio - 1:10');

% Save the optimized parameters to a .mat file
save('alfak1k3delta_opt_T10_Y3.mat', 'k1_opt', 'k3_opt', 'alfa_opt', 'delta_opt');

%%% Objective function for parameter estimation
function SSE = fun_alfak1k3delta_PS(x)

global K;

% Load experimental data
T_point = load('Data/Data_extract_041122/Time_carT10.txt');
PBT030_Y1 = load('Data/Data_extract_041122/PBT030_carT10_Y1.txt');
PBT030_Y2 = load('Data/Data_extract_041122/PBT030_carT10_Y2.txt');
PBT030_Y3 = load('Data/Data_extract_041122/PBT030_carT10_Y3.txt');

PBT030_ave = mean([PBT030_Y1, PBT030_Y2, PBT030_Y3], 2);  % Average experimental data
xdata = T_point;
ydata = PBT030_Y3;  % Use data from replicates 3 for fitting

% Assign parameters from input vector 'x'
alfa = x(1);
k1 = x(2);
k3 = x(3);
delta = x(4);
par = [alfa, K, k1, k3, delta];

% Initial conditions for the model
C_0_cellNum = ((PBT030_Y3(1) - 0.1602) / 0.0001946) / 10;
C_0 = C_0_cellNum * 0.0001946 + 0.1602;
Init = [PBT030_Y3(1); C_0];

% Solve the ODEs using ode45
[t, y] = ode45(@(t, y) CarOv(t, y, par), xdata, Init);

% Compute the sum of squared errors (SSE)
R1 = sum((ydata - y(:, 1)).^2);  % Tumor cell fit error

SSE = R1;  % Total error (SSE) for parameter fitting

end

%%% ODE function for Tumor and Car T-cell dynamics
function dydt = CarOv(t, y, par)

alfa = par(1);
K = par(2);
k1 = par(3); %theta_TC
k3 = par(4); %mu
delta = par(5);

T = y(1);  % Tumor cells
C = y(2);  % Car T-cells

% Differential equations
dTdt = alfa * T * (1 - T / K) - k1 * T * C;  % Tumor cell growth and interaction with Car T-cells
dCdt = -delta * C + k3 * T * C;  % Car T-cell decay and activation by tumor cells

dydt = [dTdt; dCdt];  % Return the derivatives

end