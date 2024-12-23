%% Phase diagram for the CAR Module

% Clear workspace and close all figures
clear all;  
close all;

% Load the required data files
load('MeanParameterValues_alfadeltak1k3.mat');  % Load average parameter values for CAR T model
load('K_opt_Tonly.mat');  % Load optimal values of K
T_point = load('Data/Data_extract_041122/Time_carT10.txt'); % Load time points (in hours)

% Set the model parameters based on the loaded data
alfa = alfa_10_ave;     % Tumor proliferation rate
K = K_opt;              % Carrying capacity K
theta1 = k1_10_ave;     % Car T killing rate (theta_TC)
theta4 = k3_10_ave;     % Car T expansion/exhaustion rate (mu)
delta = delta_10_ave;   % Car T death rate

% Compute additional parameters for the model
A = delta / alfa;  % Derived parameter A: ratio of death to proliferation rate
B = theta4 * K / alfa;  % Derived parameter B: rate involving K, mu, and alfa

% Set the time span for the ODE solver
tspan = T_point;

% Initialize concentration arrays for initial values of Tumor and CAR T cell counts
T_0 = linspace(0, 1.2, 30);  % Initial tumor concentrations from 0 to 1.2
C_0 = linspace(0, 1.2, 30);  % Initial CAR T concentrations from 0 to 1.2

% Plot the phase trajectories for different initial conditions

% Loop 1: Varying C_0, fixed T_0(end) (final value of T)
for i = 1:length(T_0)
    y0 = [T_0(end); C_0(i)];  % Set initial conditions (T_0(end) and varying C_0)
    options = odeset('reltol', 1e-8, 'abstol', [1e-8 1e-8]);  % Set ODE solver options for high accuracy
    [t, y] = ode45(@(t, y) CarOv(t, y, [A B]), tspan, y0, options);  % Solve ODEs for current initial conditions
    figure(1)
    plot(y(:, 1), y(:, 2), 'k');  % Plot the phase trajectory (T vs C)
    hold on;  % Keep the previous trajectories on the same plot
end

% Loop 2: Varying C_0, fixed T_0(2)
for i = 1:length(T_0)
    y0 = [T_0(2); C_0(i)];  % Set initial conditions (T_0(2) and varying C_0)
    options = odeset('reltol', 1e-8, 'abstol', [1e-8 1e-8]);  % Set ODE solver options
    [t, y] = ode45(@(t, y) CarOv(t, y, [A B]), tspan, y0, options);  % Solve ODEs for current initial conditions
    figure(1)
    plot(y(:, 1), y(:, 2), 'k');  % Plot the phase trajectory (T vs C)
    hold on;
end

% Loop 3: Varying T_0, fixed C_0(end)
for i = 1:length(T_0)
    y0 = [T_0(i); C_0(end)];  % Set initial conditions (varying T_0 and fixed C_0(end))
    options = odeset('reltol', 1e-8, 'abstol', [1e-8 1e-8]);  % Set ODE solver options
    [t, y] = ode45(@(t, y) CarOv(t, y, [A B]), tspan, y0, options);  % Solve ODEs for current initial conditions
    figure(1)
    plot(y(:, 1), y(:, 2), 'k');  % Plot the phase trajectory (T vs C)
    hold on;
end

% Loop 4: Varying T_0, fixed C_0(2)
for i = 1:length(T_0)
    y0 = [T_0(i); C_0(2)];  % Set initial conditions (varying T_0 and fixed C_0(2))
    options = odeset('reltol', 1e-8, 'abstol', [1e-8 1e-8]);  % Set ODE solver options
    [t, y] = ode45(@(t, y) CarOv(t, y, [A B]), tspan, y0, options);  % Solve ODEs for current initial conditions
    figure(1)
    plot(y(:, 1), y(:, 2), 'k');  % Plot the phase trajectory (T vs C)
    hold on;
end

% Define equilibrium points to mark on the plot
E1_x = T_0(2);  % Special point E1: T_0(2)
E1_y = 0;       % E1: C = 0
E2_x = 1;       % Special point E2: T = 1, C = 0
E2_y = 0;       % E2: C = 0
E3_x = T_0(2);  % Special point E3: T = A/B
E3_y = (B - A) / B;  % E3: C = (B-A)/B

% Plot equilibrium points (E1, E2, E3) with red stars
plot(E2_x, E2_y, 'r*', E3_x, E3_y, 'r*', 'MarkerSize', 10, 'LineWidth', 1.5);

% Set plot axis and labels
axis([T_0(2) 1 0 1.2]);  % Set axis limits
xlabel('T(t)', 'Interpreter', 'latex', 'FontSize', 15);  % Label x-axis (Tumor T)
ylabel('C(t)', 'Interpreter', 'latex', 'FontSize', 15);  % Label y-axis (CAR T cells C)

% Customize ticks and labels for x and y axes
xticks([T_0(2) .2 .4 .6 .8 1]);
yticks([C_0(2) .2 .4 .6 .8 1 1.2]);
xticklabels({0, 0.2, 0.4, 0.6, 0.8, 1});
yticklabels({0, 0.2, 0.4, 0.6, 0.8, 1});

% Set square axis for better visualization
axis square;

% Return at the end to stop the execution
return;

% Function definition for the ODE system (CAR T cell dynamics)
function [dydt] = CarOv(t, y, par)
    % Extract parameters A and B from the parameter vector
    A = par(1);
    B = par(2);

    % Assign state variables for Tumor (T) and CAR T cells (C)
    T = y(1);  % Tumor concentration
    C = y(2);  % CAR T cell concentration

    % Define the ODEs for the system
    dTdt = T * (1 - T) - T * C;  % Rate of change of tumor (logistic growth minus killing by CAR T cells)
    dCdt = -A * C + B * T * C;  % Rate of change of CAR T cells (death and activation by tumor)

    % Return the derivatives
    dydt = [dTdt; dCdt];
end