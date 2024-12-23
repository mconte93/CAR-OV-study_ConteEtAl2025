%% Code for plotting the bifurcation diagram of the OV module

clear all;
close all;

% Load optimized parameters from previous simulations
load('K_opt_Tonly.mat');  % Tumor carrying capacity (K)
load('Extended_Burst_size_analysis/alfabetagamma_opt_OV0008_c25_Y1.mat', 'alfa_opt', 'beta_opt', 'gamma_opt');
alfa_0008_250_Y1 = alfa_opt; beta_0008_250_Y1 = beta_opt; gamma_0008_250_Y1 = gamma_opt;
load('Extended_Burst_size_analysis/alfabetagamma_opt_OV0008_c25_Y2.mat', 'alfa_opt', 'beta_opt', 'gamma_opt');
alfa_0008_250_Y2 = alfa_opt; beta_0008_250_Y2 = beta_opt; gamma_0008_250_Y2 = gamma_opt;
load('Extended_Burst_size_analysis/alfabetagamma_opt_OV0008_c25_Y3.mat', 'alfa_opt', 'beta_opt', 'gamma_opt');
alfa_0008_250_Y3 = alfa_opt; beta_0008_250_Y3 = beta_opt; gamma_0008_250_Y3 = gamma_opt;

% Load time points for the experimental data
T_point = load('Data/Data_extract_041122/Time_OV03.txt');  % Time points (in hours)

% Calculate the average values of the parameters across the three datasets
alfa = mean([alfa_0008_250_Y1, alfa_0008_250_Y2, alfa_0008_250_Y3]);
beta = mean([beta_0008_250_Y1, beta_0008_250_Y2, beta_0008_250_Y3]);
gamma = mean([gamma_0008_250_Y1, gamma_0008_250_Y2, gamma_0008_250_Y3]);

% Constants
b = 25;  % Burst size (MOI_ICI)
omega = 0.025;  % Virus clearance rate (1/h)
K = K_opt;  % Tumor carrying capacity from previous load

% Calculate system parameters
D = gamma / alfa;
F = omega / alfa;

% Polynomial for the bifurcation diagram
p = [-1, ((F + D)^2 + D * (F + 1)), F * ((F + D)^2 + D * (F + 1)), F^3 * D];
r = roots(p);  % Solve the polynomial for equilibrium points

% Find the positive real root (Etilde)
for i = 1:length(r)
    if r(i) > 0
        Etilde = r(i);
    end
end

% Time span for the ODE integration
tspan = [0, 1000];

% Initial conditions
T_0 = 0.1;
I_0 = 0;
V_0 = 0.1;

% Array for plotting bifurcation diagram
f3 = linspace(Etilde, Etilde + 2, 30);
amp = zeros(1, length(f3));  % Store the amplitude values

% Solve ODEs for different values of E and calculate the corresponding amplitude
for i = 1:length(f3)
    y0 = [T_0; I_0; V_0];  % Initial conditions vector
    options = odeset('reltol', 1e-8, 'abstol', [1e-8, 1e-8]);
    [t, y] = ode45(@(t, y) Ov(t, y, [D, f3(i), F]), tspan, y0, options);  % Solve the system
    amp(i) = max(y(1000:end, 1));  % Take the max of the tumor population after t=1000
end

% Define the bifurcation points for plotting
f1 = [0, D * F];
f2 = [D * F, Etilde + 2];
f3 = linspace(D * F, Etilde, 30);
f30 = D * F ./ f3;  % Relation for the bifurcation diagram
f5 = linspace(Etilde, Etilde + 2, 30);
f50 = D * F ./ f5;  % Relation for the bifurcation diagram

% Plot the bifurcation diagram
figure;
plot(f1, [1, 1], 'k', f2, [1, 1], 'k--', f3, f30, 'k', f5, f50, 'k--', f3, amp, 'k.', 'LineWidth', 1, 'MarkerSize', 10);
axis([0, Etilde + 1, -0.05, 1.5]);
xticks([0, 0.28, 1, 2, 2.45, 3, 4]);
xticklabels({0, '$\bar{E}$', 1, 2, '$\tilde{E}$', 3, 4});
xlabel('$E$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$Y_3(t)$', 'Interpreter', 'latex', 'FontSize', 15);
hold on;

% Mark the equilibrium points
plot(D * F, 1, 'ro');
plot(Etilde, D * F / Etilde, 'ro');

% Add a legend
legend('Stable equilibrium', 'Unstable equilibrium', '', '', 'Stable limit cycle', 'Interpreter', 'latex', 'FontSize', 12);

% Return from the function
return;

% Define the ODE system for the CAR T cell model
function [dydt] = Ov(t, y, par)

D = par(1);  % Parameter D
E = par(2);  % Parameter E (this varies in the bifurcation diagram)
F = par(3);  % Parameter F

T = y(1);  % Tumor cells
I = y(2);  % Immune cells
V = y(3);  % Virus

% Define the ODEs
dTdt = T * (1 - T - I) - T * V;  % Tumor dynamics
dIdt = T * V - D * I;  % Immune cell dynamics
dVdt = E * I - F * V;  % Virus dynamics

dydt = [dTdt; dIdt; dVdt];  % Return the derivatives
end