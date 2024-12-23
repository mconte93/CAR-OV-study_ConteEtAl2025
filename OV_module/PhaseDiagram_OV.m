%% Phase diagrams OV module

% Initialization and Data Loading
clear all
close all

% Load the optimization parameters
load('K_opt_Tonly.mat')  % Optimized K values
load('Extended_Burst_size_analysis/alfabetagamma_opt_OV0008_c25_Y1.mat', 'alfa_opt', 'beta_opt', 'gamma_opt')  % Parameters for Y1
alfa_0008_25_Y1 = alfa_opt;
beta_0008_25_Y1 = beta_opt;
gamma_0008_25_Y1 = gamma_opt;

load('Extended_Burst_size_analysis/alfabetagamma_opt_OV0008_c25_Y2.mat', 'alfa_opt', 'beta_opt', 'gamma_opt')  % Parameters for Y2
alfa_0008_25_Y2 = alfa_opt;
beta_0008_25_Y2 = beta_opt;
gamma_0008_25_Y2 = gamma_opt;

load('Extended_Burst_size_analysis/alfabetagamma_opt_OV0008_c25_Y3.mat', 'alfa_opt', 'beta_opt', 'gamma_opt')  % Parameters for Y3
alfa_0008_25_Y3 = alfa_opt;
beta_0008_25_Y3 = beta_opt;
gamma_0008_25_Y3 = gamma_opt;

% Averaging the parameters across all three conditions (Y1, Y2, Y3)
alfa = mean([alfa_0008_25_Y1, alfa_0008_25_Y2, alfa_0008_25_Y3]);
beta = mean([beta_0008_25_Y1, beta_0008_25_Y2, beta_0008_25_Y3]);
gamma = mean([gamma_0008_25_Y1, gamma_0008_25_Y2, gamma_0008_25_Y3]);

% Additional parameters
b = 25;  % Constant for the virus dynamics
omega = 0.025;  % Constant for the virus decay rate
K = K_opt;  % Load K optimization

% Load time data for the system
T_point = load('Data/Data_extract_041122/Time_OV03.txt');  % Time points in hours

%% Derived Parameters
D = gamma / alfa;  % Derived parameter D
E = 3;  % E parameter: set to 3 for this case
F = omega / alfa;  % Derived parameter F

% Output the derived parameters to verify
disp(D * F)
disp(E)

%% Function \mathcal{H} for further analysis
% Function to derive the Hcal curve
E_vec = linspace(0, 10, 100);  % Vector for E values

% Define the mathematical expression for Hcal
Hcal = -E_vec.^3 + E_vec.^2 .* ((F + D).^2 + D .* (F + 1)) + F .* E_vec .* ((F + D).^2 + D .* (F + 1)) + F^3 .* D;

% Plot Hcal as a function of E_vec
figure(2)
plot(E_vec, Hcal, [0 10], [0 0], 'r')  % Plot Hcal curve with reference line at y=0

%% Simulation with ODEs (Trajectory Plots in Phase Space)
tspan = T_point;  % Define time span

% Initial conditions for various cases (T_0, I_0, V_0)
T_0 = linspace(0, 1, 5);
I_0 = linspace(0, 1, 5);
V_0 = linspace(0, 1, 5);

% Loop through initial conditions and solve ODEs for each case
for i = 1:length(I_0)
    for j = 1:length(V_0)
        % Initial state vector
        y0 = [T_0(end); I_0(i); V_0(j)];
        
        % Solve ODE using ode45 with specified options
        options = odeset('reltol', 1e-8, 'abstol', [1e-8, 1e-8]);
        [t, y] = ode45(@(t, y) Ov(t, y, [D, E, F]), [], y0, options);
        
        % Plot phase space trajectory (T vs I vs V)
        figure(1)
        plot3(y(:, 1), y(:, 2), y(:, 3), 'k');  % Phase space trajectory in 3D
        hold on
    end
end

% Repeat the process for other initial conditions
for i = 1:length(I_0)
    for j = 1:length(V_0)
        y0 = [T_0(1); I_0(i); V_0(j)];
        options = odeset('reltol', 1e-8, 'abstol', [1e-8, 1e-8]);
        [t, y] = ode45(@(t, y) Ov(t, y, [D, E, F]), tspan, y0, options);
        figure(1)
        plot3(y(:, 1), y(:, 2), y(:, 3), 'k');  % Phase space trajectory
        hold on
    end
end

% Further simulations with different initial conditions for T_0, I_0, V_0
for i = 1:length(T_0)
    for j = 1:length(I_0)
        y0 = [T_0(i); I_0(j); V_0(end)];
        options = odeset('reltol', 1e-8, 'abstol', [1e-8, 1e-8]);
        [t, y] = ode45(@(t, y) Ov(t, y, [D, E, F]), tspan, y0, options);
        figure(1)
        plot3(y(:, 1), y(:, 2), y(:, 3), 'k');  % Phase space trajectory
        hold on
    end
end

% More simulations with initial conditions for V_0(1) and V_0(end)
for i = 1:length(T_0)
    for j = 1:length(I_0)
        y0 = [T_0(i); I_0(end); V_0(j)];
        options = odeset('reltol', 1e-8, 'abstol', [1e-8, 1e-8]);
        [t, y] = ode45(@(t, y) Ov(t, y, [D, E, F]), tspan, y0, options);
        figure(1)
        plot3(y(30:end, 1), y(30:end, 2), y(30:end, 3), 'k');  % Plot phase trajectory from t=30 onwards
        hold on
    end
end

%% Mark Specific Points in Phase Space
% Marking equilibrium points in the phase space trajectory
E1_x = 0;
E1_y = 0;
E1_z = 0;
E2_x = 1;
E2_y = 0;
E2_z = 0;
E3_x = D * F / E;
E3_y = F / E * (E - D * F) / (E + F);
E3_z = (E - D * F) / (E + F);

% Plot these points on the trajectory plot
figure(1)
plot3(E1_x, E1_y, E1_z, 'r*', E2_x, E2_y, E2_z, 'r*', E3_x, E3_y, E3_z, 'r*', 'MarkerSize', 10, 'LineWidth', 1.5)

% Customize the plot
xlabel('$Y_1(t)$', 'Interpreter', 'latex', 'FontSize', 15)
ylabel('$Y_2(t)$', 'Interpreter', 'latex', 'FontSize', 15)
zlabel('$Y_3(t)$', 'Interpreter', 'latex', 'FontSize', 15)
axis square

%% ODE System Definition
% Define the system of ODEs for the model
function [dydt] = Ov(t, y, par)
    % Parameters for the system
    D = par(1);  % Parameter D
    E = par(2);  % Parameter E
    F = par(3);  % Parameter F

    % State variables
    T = y(1);  % Tumor cells
    I = y(2);  % Infected cells
    V = y(3);  % Virus concentration

    % Define the differential equations for the system
    dTdt = T * (1 - T - I) - T * V;  % Tumor growth and interaction with virus
    dIdt = T * V - D * I;  % Infected cell dynamics
    dVdt = E * I - F * V;  % Virus dynamics

    % Return the derivatives as a column vector
    dydt = [dTdt; dIdt; dVdt];
end