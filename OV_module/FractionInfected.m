%% Code for infected fraction analysis
clear all
close all

% Load the Fraction Data (contains data for fractions of infected cells)
load('Fraction_Data.mat')

% Figure for fractions of infected cells at two different times
figure(3)

%% Plot 1: Fractions at T=24 hours (MOI data)
subplot(1,2,1)
semilogx(x, Frac24_ave, 'r-o', 'MarkerSize', 6)  % Plot average fraction at T=24h
hold on
errorbar(x, Frac24_ave, Std_frac24, Std_frac24, 'r-', 'LineWidth', 1)  % Error bars for T=24h

% Interpolated data points for specific MOIs
semilogx(0.0008, Int1(2), 'r*', 0.002, Int1(4), 'r*', 0.03, Int1(6), 'r*', 'MarkerSize', 8)

% Set plot labels, title, and axis limits
title('T=24h', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('MOI', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('Fractions of infected cells', 'Interpreter', 'latex', 'FontSize', 16)
ylim([0 1.5])
xlim([0.0005 1])
axis square

%% Plot 2: Fractions at T=48 hours (MOI data)
subplot(1,2,2)
semilogx(x, Frac48_ave, 'b-o', 'MarkerSize', 6)  % Plot average fraction at T=48h
hold on
errorbar(x, Frac48_ave, Std_frac48, Std_frac48, 'b-', 'LineWidth', 1)  % Error bars for T=48h

% Interpolated data points for specific MOIs
semilogx(0.0008, Int2(2), 'b*', 0.002, Int2(4), 'b*', 0.03, Int2(6), 'b*', 'MarkerSize', 8)

% Set plot labels, title, and axis limits
legend('Data', '', 'Interp. data', 'Interpreter', 'latex', 'FontSize', 16)
title('T=48h', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('MOI', 'Interpreter', 'latex', 'FontSize', 16)
ylim([0 1.5])
xlim([0.0005 1])
axis square

%% Plot 3: Ratio of Fractions at T=48h / T=24h
figure(4)
semilogx(x, Ratio, 'k-o', 'MarkerSize', 6, 'LineWidth', 1)  % Plot ratio of T=48h to T=24h
hold on
errorbar(x, Ratio, Std_ratio, Std_ratio, 'k-', 'LineWidth', 0.6)  % Error bars for ratio

% Interpolated data points for specific MOIs
semilogx(0.0008, Int2(2) / Int1(2), 'k*', 0.002, Int2(4) / Int1(4), 'k*', 0.03, Int2(6) / Int1(6), 'k*', 'MarkerSize', 8)

% Set plot labels, title, and axis limits
ylabel('$\frac{I(T=48)}{I(T=24)}$', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('MOI', 'Interpreter', 'latex', 'FontSize', 16)
legend('Data', '', 'Interp. Data', 'Interpreter', 'latex', 'FontSize', 16)
title('Ratio', 'Interpreter', 'latex', 'FontSize', 16)
xlim([0.0005-0.00001 1])
axis square

%% Forward Model for MOI = 0.03 (Numerical Results)

% Load data for MOI=0.03 and relevant parameters
T_point = load('Data/Data_extract_041122/Time_OV03.txt');  % Time data (hours)
load ('K_opt_Tonly.mat', 'K_opt')  % Load optimized K values
load('Extended_Burst_size_analysis/alfabetagamma_opt_OV03_c25.mat');  % Load optimal parameters for MOI=0.03

% Load experimental data for MOI=0.03
PBT030_Y1 = load('Data_extract_041122/PBT030_OV03_Y1.txt');  % Data at T=0
PBT030_Y2 = load('Data/Data_extract_041122/PBT030_OV03_Y2.txt');  % Data at T=24
PBT030_Y3 = load('Data/Data_extract_041122/PBT030_OV03_Y3.txt');  % Data at T=48

% Average the experimental data for MOI=0.03
PBT030_ave = mean([PBT030_Y1 PBT030_Y2 PBT030_Y3], 2);

% Set parameters for the model
K = K_opt;
alfa = alfa_opt;
gamma = gamma_opt;
beta = beta_opt;
omega = 0.5;
b = 25;

% Initial conditions (initial values for T, I, V)
time = T_point;
T_0 = PBT030_ave(1);
I_0 = 0;
V_0_03 = 0.03;
Init = [T_0; I_0; V_0_03];

% Solve the ODEs using ode45
[t, y] = ode45(@(t, y) Ov(t, y, [alfa K beta gamma b omega]), time, Init, []);

% Find the indices for T=24h and T=48h
indx_24 = find(fix(t) == 24);
indx_24 = indx_24(1);
indx_48 = find(fix(t) == 48);
indx_48 = indx_48(1);

% Calculate the infected fraction at T=24h and T=48h
I24_03 = y(indx_24, 2) / (y(indx_24, 1) + y(indx_24, 2));
I48_03 = y(indx_48, 2) / (y(indx_48, 1) + y(indx_48, 2));

% Plot the model results for MOI = 0.03
figure(3)
subplot(1, 2, 1)
hold on
semilogx(V_0_03, I24_03, 'k*', 'MarkerSize', 8)
subplot(1, 2, 2)
semilogx(V_0_03, I48_03, 'k*', 'MarkerSize', 8)

% Plot the ratio of infected fractions
figure(4)
semilogx(V_0_03, I48_03 / I24_03, 'r*', 'MarkerSize', 6)
hold on

%% Forward Model for MOI = 0.002 (Similar to MOI = 0.03)

% Load the data for MOI=0.002 and other relevant parameters
T_point = load('Data/Data_extract_041122/Time_OV002.txt');
load('Extended_Burst_size_analysis/alfabetagamma_opt_OV002_c25.mat');  % MOI = 0.002 parameters
PBT030_Y1 = load('Data/Data_extract_041122/PBT030_OV002_Y1.txt');
PBT030_Y2 = load('Data/Data_extract_041122/PBT030_OV002_Y2.txt');
PBT030_Y3 = load('Data/Data_extract_041122/PBT030_OV002_Y3.txt');
PBT030_ave = mean([PBT030_Y1 PBT030_Y2 PBT030_Y3], 2);

% Set parameters and initial conditions
V_0_002 = 0.002;
Init = [PBT030_ave(1); 0; V_0_002];  % Initial conditions for MOI=0.002
[t, y] = ode45(@(t, y) Ov(t, y, [alfa K beta gamma b omega]), time, Init, []);

% Calculate infected fractions at T=24h and T=48h
indx_24 = find(fix(t) == 24);
indx_48 = find(fix(t) == 48);
I24_002 = y(indx_24, 2) / (y(indx_24, 1) + y(indx_24, 2));
I48_002 = y(indx_48, 2) / (y(indx_48, 1) + y(indx_48, 2));

% Plot the model results for MOI = 0.002
figure(3)
subplot(1, 2, 1)
hold on
semilogx(V_0_002, I24_002, 'k*', 'MarkerSize', 8)
subplot(1, 2, 2)
semilogx(V_0_002, I48_002, 'k*', 'MarkerSize', 8)

% Plot the ratio of infected fractions
figure(4)
semilogx(V_0_002, I48_002 / I24_002, 'r*', 'MarkerSize', 6)
hold on

%% Forward Model for MOI = 0.0008 (Similar to MOI = 0.03)

% Load the data for MOI=0.0008 and other relevant parameters
T_point = load('Data/Data_extract_041122/Time_OV0008.txt');
load('Extended_Burst_size_analysis/alfabetagamma_opt_OV0008_c25.mat');
PBT030_Y1 = load('Data/Data_extract_041122/PBT030_OV0008_Y1.txt');
PBT030_Y2 = load('Data/Data_extract_041122/PBT030_OV0008_Y2.txt');
PBT030_Y3 = load('Data/Data_extract_041122/PBT030_OV0008_Y3.txt');
PBT030_ave = mean([PBT030_Y1 PBT030_Y2 PBT030_Y3], 2);

% Set parameters and initial conditions for MOI = 0.0008
V_0_0008 = 0.0008;
Init = [PBT030_ave(1); 0; V_0_0008];
[t, y] = ode45(@(t, y) Ov(t, y, [alfa K beta gamma b omega]), time, Init, []);

% Calculate infected fractions at T=24h and T=48h
indx_24 = find(fix(t) == 24);
indx_48 = find(fix(t) == 48);
I24_0008 = y(indx_24, 2) / (y(indx_24, 1) + y(indx_24, 2));
I48_0008 = y(indx_48, 2) / (y(indx_48, 1) + y(indx_48, 2));

% Plot the model results for MOI = 0.0008
figure(3)
subplot(1, 2, 1)
hold on
semilogx(V_0_0008, I24_0008, 'k*', 'MarkerSize', 8)
legend('Data', '', 'Interp. data', '', '', 'Num. Res.', '', '', 'Interpreter', 'latex', 'FontSize', 16)
subplot(1, 2, 2)
semilogx(V_0_0008, I48_0008, 'k*', 'MarkerSize', 8)
legend('Data', '', 'Interp. data', '', '', 'Num. Res.', '', '', 'Interpreter', 'latex', 'FontSize', 16)

% Plot the ratio of infected fractions
figure(4)
semilogx(V_0_0008, I48_0008 / I24_0008, 'r*', 'MarkerSize', 6)
hold on

%% Save Results
save(['NumRes_Inf_b_', num2str(b), '_omega_', num2str(omega), '.mat'], ...
     'V_0_03', 'I24_03', 'I48_03', 'V_0_002', 'I24_002', 'I48_002', 'V_0_0008', 'I24_0008', 'I48_0008')

%% ODE Model Function
function [dydt] = Ov(t, y, par)
    % Define parameters
    alfa = par(1);
    K = par(2);
    beta = par(3);
    gamma = par(4);
    b = par(5);
    omega = par(6);
    
    % State variables
    T = y(1);  % Tumor cells
    I = y(2);  % Infected cells
    V = y(3);  % Virus concentration
    
    % Define the differential equations
    dTdt = alfa * T * (1 - (T + I) / K) - beta * T * V;  % Tumor growth
    dIdt = beta * T * V - gamma * I;  % Infected cell dynamics
    dVdt = b * gamma * I - omega * V;  % Virus dynamics
    
    % Return the derivatives
    dydt = [dTdt; dIdt; dVdt];
end