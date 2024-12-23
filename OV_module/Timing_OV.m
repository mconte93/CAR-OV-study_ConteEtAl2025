%% Study of therapy timing administration for the OV monotherapy case

%% Initialization

clear all
% close all

global b

% Load time points and data for the treatment groups (in hours)
T_point = load('Data/Data_extract_041122/Time_OV03.txt');  % Time points in hours
PBT030_Y1 = load('Data/Data_extract_041122/PBT030_OV03_Y1.txt');  % Data for treatment group 1 (hours)
PBT030_Y2 = load('Data/Data_extract_041122/PBT030_OV03_Y2.txt');  % Data for treatment group 2 (hours)
PBT030_Y3 = load('Data/Data_extract_041122/PBT030_OV03_Y3.txt');  % Data for treatment group 3 (hours)

% Average of the three treatment groups
PBT030_ave = mean([PBT030_Y1 PBT030_Y2 PBT030_Y3], 2);

% Load model parameters (from pre-computed data)
load('K_opt_Tonly.mat', 'K_opt')
K = K_opt;  % Carry capacity (K)

% Load infection-related parameters (from pre-computed data)
load('Extended_Burst_size_analysis/alfabetagamma_opt_OV03_c25.mat', 'gamma_opt', 'beta_opt', 'alfa_opt')
alfa = alfa_opt;  % Growth rate of tumor cells
beta = beta_opt;  % Infection rate (viral spread)
gamma = gamma_opt;  % Recovery rate of infected cells
omega = 0.05;  % Virus clearance rate
b = 25;  % Burst size (100 PFU/cell)

% Define parameter vector for ease of passing into ODE function
par = [alfa, K, beta, gamma, b, omega];

%% Main loop: Solve ODE for each time point (0 to 38)

% Initialize arrays to store results
Area_OV = zeros(1, 39);
FinalT_OV = zeros(1, 39);
Y = zeros(1, 39);

% Loop through time points (k=0 to 38)
for k = 0:38
    Y(k + 1) = T_point(1 + k * 4);  % Viral administration time
    adm_time_V = 1 + k * 4;  % Corresponding index for the administration time
    
    if adm_time_V > 1
        % Time before viral administration
        time_1 = T_point(1:adm_time_V - 1);
        T_0 = PBT030_ave(1);  % Initial tumor size
        I_0 = 0;  % Initial infected cells
        V_0 = 0;  % Initial viral load
        Init_1 = [T_0; I_0; V_0];  % Initial conditions for ODE
        
        % Solve ODE for pre-administration period
        [t_1, y_1] = ode45(@(t, y) Ov(t, y, par), time_1, Init_1, []);
        
        % Time after viral administration
        time_2 = T_point(adm_time_V:end);
        T_0 = y_1(end, 1);  % Initial tumor size from previous solution
        I_0 = 0;  % No infected cells at the beginning of this phase
        V_0 = 0.03;  % Initial viral load after administration
        Init_2 = [T_0; I_0; V_0];  % Initial conditions for ODE after viral admin
        
        % Solve ODE for post-administration period
        [t_2, y_2] = ode45(@(t, y) Ov(t, y, par), time_2, Init_2, []);
        
        % Calculate the area under the curve (numerical integration)
        A_num = 0;
        t_area = [t_1; t_2(2:end)];  % Combine time vectors from both phases
        y_area = [y_1(:, 1:2); y_2(2:end, 1:2)];  % Combine tumor size and infected cells
        
        % Trapezoidal rule for numerical integration
        for s = 1:(size(t_area, 1) - 1)
            h = (t_area(s + 1) - t_area(s));  % Time step
            b1 = y_area(s, 1) + y_area(s, 2);  % Sum of tumor and infected cells at time t_s
            b2 = y_area(s + 1, 1) + y_area(s + 1, 2);  % Sum of tumor and infected cells at time t_s+1
            A_num = A_num + (b1 + b2) * h / 2;  % Trapezoidal rule integration
        end
        
        Area_OV(k + 1) = A_num;  % Store area under curve
        FinalT_OV(k + 1) = y_area(end, 1) + y_area(end, 2);  % Final tumor size (sum of tumor and infected cells)
    else
        % No viral administration (use the entire time series)
        time_1 = T_point(1:end);
        T_0 = PBT030_ave(1);
        I_0 = 0;
        V_0 = 0.03;
        Init_1 = [T_0; I_0; V_0];  % Initial conditions for ODE
        
        % Solve ODE for the entire time series (no administration phase)
        [t_1, y_1] = ode45(@(t, y) Ov(t, y, par), time_1, Init_1, []);
        
        % Calculate the area under the curve
        A_num = 0;
        t_area = t_1;
        y_area = y_1(:, 1:2);
        
        % Trapezoidal rule for numerical integration
        for s = 1:(size(t_area, 1) - 1)
            h = (t_area(s + 1) - t_area(s));  % Time step
            b1 = y_area(s, 1) + y_area(s, 2);  % Sum of tumor and infected cells at time t_s
            b2 = y_area(s + 1, 1) + y_area(s + 1, 2);  % Sum of tumor and infected cells at time t_s+1
            A_num = A_num + (b1 + b2) * h / 2;  % Trapezoidal rule integration
        end
        
        Area_OV(k + 1) = A_num;  % Store area under curve
        FinalT_OV(k + 1) = y_area(end, 1) + y_area(end, 2);  % Final tumor size
    end
end

%% Plotting the results

% Create figure to plot results
figure(1)

% Plot AUC vs viral administration time
subplot(1, 2, 1)
plot(Y, Area_OV(1, :), 'k-.')
hold on
xlim([min(Y) max(Y)])
xlabel('Viral administration time [h]', 'Interpreter', 'latex', 'FontSize', 13)
ylabel('AUC', 'Interpreter', 'latex', 'FontSize', 13)
axis square

% Plot final tumor value vs viral administration time
subplot(1, 2, 2)
plot(Y, FinalT_OV(1, :), 'k-.')
hold on
xlim([min(Y) max(Y)])
xlabel('Viral administration time [h]', 'Interpreter', 'latex', 'FontSize', 13)
ylabel('Final tumor value', 'Interpreter', 'latex', 'FontSize', 13)
axis square

%% Save results
save('Timing_OV03.mat', 'Y', 'Area_OV', 'FinalT_OV')


%% ODE function that defines the system of differential equations

function [dydt] = Ov(t, y, par)
    % Unpack the parameter vector
    alfa = par(1);
    K = par(2);
    beta = par(3);
    gamma = par(4);
    b = par(5);
    omega = par(6);

    % Define the state variables
    T = y(1);  % Tumor size
    I = y(2);  % Infected cells
    V = y(3);  % Virus load

    % Differential equations
    dTdt = alfa * T * (1 - (T + I) / K) - beta * T * V;  % Tumor growth and virus-induced death
    dIdt = beta * T * V - gamma * I;  % Infected cells dynamics
    dVdt = b * gamma * I - omega * V;  % Virus dynamics

    % Return the derivatives
    dydt = [dTdt; dIdt; dVdt];
end