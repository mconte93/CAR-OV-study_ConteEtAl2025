%% Study on therapy administration for monotherapy with CAR T cells only
clear all;
close all;

% Load the data files for time points and CAR T cell dynamics
T_point = load('Data/Data_extract_041122/Time_carT10.txt');  % Time points (in hours)
PBT030_Y1 = load('Data/Data_extract_041122/PBT030_carT10_Y1.txt');  % Data for PBT030 (Y1)
PBT030_Y2 = load('Data/Data_extract_041122/PBT030_carT10_Y2.txt');  % Data for PBT030 (Y2)
PBT030_Y3 = load('Data/Data_extract_041122/PBT030_carT10_Y3.txt');  % Data for PBT030 (Y3)

% Calculate average data from the three datasets
PBT030_ave = mean([PBT030_Y1, PBT030_Y2, PBT030_Y3], 2);  % Averaging the three data series

% Load model parameters from external files
load('K_opt_Tonly.mat', 'K_opt');  % Optimal value of K (Tumor carrying capacity)
K = K_opt;  % Assign K

% Load optimized parameters for the model
load('alfak1k3delta_opt_T10.mat', 'delta_opt', 'k3_opt', 'k1_opt', 'alfa_opt');
delta = delta_opt;  % Tumor cell death rate
theta4 = k3_opt;    % Tumor expansion/exhaustion rate
alfa = alfa_opt;    % Tumor growth rate
theta1 = k1_opt;    % Tumor-killing rate by CAR T cells

% Pack parameters into a vector for passing to the ODE function
par = [alfa, K, theta1, theta4, delta];

% Loop over different CAR T administration times (from 0 to 38)
for i = 0:38
    % Define the time points for each administration
    X(i+1) = T_point(1 + i * 4);  % T_point is assumed to be sampled every 4 hours
    adm_time_C = 1 + i * 4;  % Administration time index
    
    % If administration time is greater than 1, break the time into two parts
    if adm_time_C > 1
        % First part: before CAR T administration
        time_1 = T_point(1:adm_time_C - 1);  % Time points before the administration
        T_0 = PBT030_ave(1);  % Initial tumor concentration (first value in the average data)
        C_0 = 0;  % Initial CAR T cell concentration (zero before administration)
        Init_1 = [T_0; C_0];  % Initial conditions
        
        % Solve the ODEs for the first part (before CAR T)
        [t_1, y_1] = ode45(@(t, y) Car(t, y, par), time_1, Init_1, []);

        % Second part: after CAR T administration
        time_2 = T_point(adm_time_C - 1:end);  % Time points after administration
        T_0 = PBT030_ave(1);  % Use the same initial tumor concentration
        % Compute the initial CAR T cell concentration based on the tumor level
        C_0_cellNum = ((T_0 - 0.1602) / 0.0001946) / 10;
        C_0 = C_0_cellNum * 0.0001946 + 0.1602;  % Adjusted CAR T cell concentration
        T_0 = y_1(end, 1);  % Use the last tumor value from the first simulation
        Init_2 = [T_0; C_0];  % Initial conditions for the second part
        
        % Solve the ODEs for the second part (after CAR T administration)
        [t_2, y_2] = ode45(@(t, y) Car(t, y, par), time_2, Init_2, []);
        
        % Calculate the total area under the curve (AUC) for the tumor
        A_num = 0;
        t_area = [t_1; t_2(2:end)];  % Concatenate the time arrays (excluding the first point of the second part)
        y_area = [y_1(:, 1:2); y_2(2:end, 1:2)];  % Concatenate the solution arrays (excluding the first point of the second part)

        % Calculate the AUC using the trapezoidal rule
        for s = 1:(size(t_area, 1) - 1)
            h = (t_area(s + 1) - t_area(s));  % Time step
            b1 = y_area(s, 1);  % Tumor value at time t(s)
            b2 = y_area(s + 1, 1);  % Tumor value at time t(s+1)
            A_num = A_num + (b1 + b2) * h / 2;  % Trapezoidal area calculation
        end
        
        % Store the calculated area and final tumor value
        Area_CAR(i + 1) = A_num;
        FinalT_CAR(i + 1) = y_area(end, 1);
    else
        % Case when administration time is 1, only one time period
        time_1 = T_point(1:end);  % Use the entire time span
        T_0 = PBT030_ave(1);  % Initial tumor concentration
        % Compute the initial CAR T cell concentration as before
        C_0_cellNum = ((T_0 - 0.1602) / 0.0001946) / 10;
        C_0 = C_0_cellNum * 0.0001946 + 0.1602;  % Adjusted CAR T cell concentration
        Init_1 = [T_0; C_0];  % Initial conditions

        % Solve the ODEs for the entire time span
        [t_1, y_1] = ode45(@(t, y) Car(t, y, par), time_1, Init_1, []);
        
        % Calculate the total area under the curve (AUC) for the tumor
        A_num = 0;
        t_area = t_1;
        y_area = y_1(:, 1);
        
        % Calculate the AUC using the trapezoidal rule
        for s = 1:(size(t_area, 1) - 1)
            h = (t_area(s + 1) - t_area(s));  % Time step
            b1 = y_area(s);  % Tumor value at time t(s)
            b2 = y_area(s + 1);  % Tumor value at time t(s+1)
            A_num = A_num + (b1 + b2) * h / 2;  % Trapezoidal area calculation
        end
        
        % Store the calculated area and final tumor value
        Area_CAR(i + 1) = A_num;
        FinalT_CAR(i + 1) = y_area(end);
    end
end

% Plot the results in two subplots: AUC vs Administration time, and Final Tumor Value vs Administration time
figure(1)
subplot(1, 2, 1)
plot(X, Area_CAR, 'k-.')
xlim([min(X), max(X)])
xlabel('CAR T administration [h]', 'Interpreter', 'latex', 'FontSize', 13)
ylabel('Total tumor area', 'Interpreter', 'latex', 'FontSize', 13)
axis square
hold on

subplot(1, 2, 2)
plot(X, FinalT_CAR, 'k-.')
xlim([min(X), max(X)])
xlabel('CAR T administration [h]', 'Interpreter', 'latex', 'FontSize', 13)
ylabel('Final tumor value', 'Interpreter', 'latex', 'FontSize', 13)
axis square
hold on

% Save the results to a file
save('Timing_CAR10.mat', 'X', 'Area_CAR', 'FinalT_CAR')

% Define the Car function (ODE system) that models the tumor and CAR T cell dynamics
function [dydt] = Car(t, y, par)
    % Extract parameters
    alfa = par(1);
    K = par(2);
    k1 = par(3);
    k3 = par(4);
    delta = par(5);
    
    % Define state variables: Tumor (T) and CAR T cells (C)
    T = y(1);  % Tumor concentration
    C = y(2);  % CAR T concentration
    
    % Define the ODEs
    dTdt = alfa * T * (1 - (T / K)) - k1 * T * C;  % Tumor dynamics (logistic growth minus CAR T killing)
    dCdt = -delta * C + k3 * T * C;  % CAR T cell dynamics (death and activation by tumor)

    % Return the derivatives
    dydt = [dTdt; dCdt];
end
