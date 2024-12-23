%% Combined therapy timing administration analysis

% Clear workspace and close all figures
clear all;
close all;

% Global variables
global b;

% Load experimental data
T_point = load('Data/Data_extract_041122/Time_50_002.txt');      % Time points in hours
PBT030_Y1 = load('Data/Data_extract_041122/PBT030_50_002_Y1.txt'); % Experimental data Y1 (tumor size)
PBT030_Y2 = load('Data/Data_extract_041122/PBT030_50_002_Y2.txt'); % Experimental data Y2 (tumor size)
PBT030_Y3 = load('Data/Data_extract_041122/PBT030_50_002_Y3.txt'); % Experimental data Y3 (tumor size)
PBT030_ave = mean([PBT030_Y1, PBT030_Y2, PBT030_Y3], 2);         % Average of tumor size data

% Load optimized parameters for the tumor model
load('K_opt_Tonly.mat', 'K_opt');  % Carrying capacity K
K = K_opt;

% Load optimized parameters for delta, k3, and k1
load('alfak1k3delta_opt_T50.mat', 'delta_opt', 'k3_opt', 'k1_opt');
delta = delta_opt; 
theta4 = k3_opt;

% Load optimized parameters for gamma and beta
load('alfabetagamma_opt_OV002_c25.mat', 'gamma_opt', 'beta_opt');
gamma = gamma_opt; 
omega = 0.05;  % Viral clearance rate
b = 25;         % Burst size

% Load average values for CAR T-cell and tumor model parameters
load('MeanParameterValues_theta123alfabeta.mat', ...
    'beta_50_002_ave', 'theta1_50_002_ave', 'alfa_50_002_ave', ...
    'theta2_50_002_ave', 'theta3_50_002_ave');

beta = beta_50_002_ave;
alfa = alfa_50_002_ave;
theta1 = theta1_50_002_ave;
theta2 = theta2_50_002_ave;
theta3 = theta3_50_002_ave;

% Combine all parameters into a single vector for easier use in ODE
par = [alfa, K, beta, theta1, gamma, theta2, b, theta3, omega, delta, theta4];

% Initialize matrices to store results
Area_comb = zeros(39, 39);   % Total tumor area
FinalT_comb = zeros(39, 39); % Final tumor size

% Simulation loop: V and CAR T-cell administration times
for i=0:38
X(i+1)=T_point(1+i*4);
adm_time_C=1+i*4;  

for k=0:38
Y(k+1)=T_point(1+k*4);
adm_time_V=1+k*4;

if adm_time_C>adm_time_V % Case 1: CAR T-cell admin after Virus admin

if adm_time_V>1
% First phase: Before virus administration
time_1=T_point(1:adm_time_V-1);
T_0=PBT030_ave(1);
C_0=0;
I_0=0;
V_0=0;
Init_1=[T_0;I_0;V_0;C_0];

[t_1,y_1] = ode45(@(t,y) TumorCarOv(t,y,par),time_1,Init_1,[]);

% Second phase: After virus administration, before CAR T-cell
time_2=T_point(adm_time_V-1:adm_time_C-1);
T_0=y_1(end,1);
C_0=0;
I_0=0;
V_0=0.002;
Init_2=[T_0;I_0;V_0;C_0];

[t_2,y_2] = ode45(@(t,y) TumorCarOv(t,y,par),time_2,Init_2,[]);

% Third phase: After CAR T-cell administration
time_3=T_point(adm_time_C-1:end);
T_0=PBT030_ave(1);
C_0_cellNum=((T_0-0.1602)/0.0001946)/25;
C_0=C_0_cellNum*0.0001946+0.1602;
T_0=y_2(end,1);
I_0=y_2(end,2);
V_0=y_2(end,3);
Init_3=[T_0;I_0;V_0;C_0];

[t_3,y_3] = ode45(@(t,y) TumorCarOv(t,y,par),time_3,Init_3,[]);

% Calculate AUC for tumor dynamics (area integration)
A_num=0;
t_area=[t_1;t_2(2:end);t_3(2:end)];
y_area=[y_1(:,1:2);y_2(2:end,1:2);y_3(2:end,1:2)];
for s=1:(size(t_area,1)-1)
    h=(t_area(s+1)-t_area(s));
    b1=y_area(s,1)+y_area(s,2);
    b2=y_area(s+1,1)+y_area(s+1,2);
    A_num=A_num+(b1+b2)*h/2;
end
Area_comb(i+1,k+1)=A_num;
FinalT_comb(i+1,k+1)=y_area(end,1)+y_area(end,2);

    else
time_1=T_point(1:adm_time_C-1);
T_0=PBT030_ave(1);
C_0=0;
I_0=0;
V_0=0.002;
Init_1=[T_0;I_0;V_0;C_0];

[t_1,y_1] = ode45(@(t,y) TumorCarOv(t,y,par),time_1,Init_1,[]);

time_2=T_point(adm_time_C-1:end);
T_0=PBT030_ave(1);
C_0_cellNum=((T_0-0.1602)/0.0001946)/25;
C_0=C_0_cellNum*0.0001946+0.1602;
T_0=y_1(end,1);
I_0=y_1(end,2);
V_0=y_1(end,3);
Init_2=[T_0;I_0;V_0;C_0];

[t_2,y_2] = ode45(@(t,y) TumorCarOv(t,y,par),time_2,Init_2,[]);

A_num=0;
t_area=[t_1;t_2(2:end)];
y_area=[y_1(:,1:2);y_2(2:end,1:2)];
for s=1:(size(t_area,1)-1)
    h=(t_area(s+1)-t_area(s));
    b1=y_area(s,1)+y_area(s,2);
    b2=y_area(s+1,1)+y_area(s+1,2);
    A_num=A_num+(b1+b2)*h/2;
end
Area_comb(i+1,k+1)=A_num;
FinalT_comb(i+1,k+1)=y_area(end,1)+y_area(end,2);
end

elseif adm_time_C<adm_time_V  % Similar structure but with CAR T administration first (CAR T-cell only)

    if adm_time_C>1
time_1=T_point(1:adm_time_C-1);
T_0=PBT030_ave(1);
C_0=0;
I_0=0;
V_0=0;
Init_1=[T_0;I_0;V_0;C_0];

[t_1,y_1] = ode45(@(t,y) TumorCarOv(t,y,par),time_1,Init_1,[]);

time_2=T_point(adm_time_C-1:adm_time_V-1);
T_0=PBT030_ave(1);
C_0_cellNum=((T_0-0.1602)/0.0001946)/25;
C_0=C_0_cellNum*0.0001946+0.1602;
T_0=y_1(end,1);
I_0=0;
V_0=0;
Init_2=[T_0;I_0;V_0;C_0];

[t_2,y_2] = ode45(@(t,y) TumorCarOv(t,y,par),time_2,Init_2,[]);

time_3=T_point(adm_time_V-1:end);
T_0=y_2(end,1);
I_0=0;
V_0=0.002;
C_0=y_2(end,4);
Init_3=[T_0;I_0;V_0;C_0];

[t_3,y_3] = ode45(@(t,y) TumorCarOv(t,y,par),time_3,Init_3,[]);

A_num=0;
t_area=[t_1;t_2(2:end);t_3(2:end)];
y_area=[y_1(:,1:2);y_2(2:end,1:2);y_3(2:end,1:2)];
for s=1:(size(t_area,1)-1)
    h=(t_area(s+1)-t_area(s));
    b1=y_area(s,1)+y_area(s,2);
    b2=y_area(s+1,1)+y_area(s+1,2);
    A_num=A_num+(b1+b2)*h/2;
end
Area_comb(i+1,k+1)=A_num;
FinalT_comb(i+1,k+1)=y_area(end,1)+y_area(end,2);

else 
       
time_1=T_point(1:adm_time_V-1);
T_0=PBT030_ave(1);
C_0_cellNum=((T_0-0.1602)/0.0001946)/25;
C_0=C_0_cellNum*0.0001946+0.1602;
I_0=0;
V_0=0;
Init_1=[T_0;I_0;V_0;C_0];

[t_1,y_1] = ode45(@(t,y) TumorCarOv(t,y,par),time_1,Init_1,[]);

time_2=T_point(adm_time_V-1:end);
T_0=y_1(end,1);
I_0=0;
V_0=0.002;
C_0=y_1(end,4);
Init_2=[T_0;I_0;V_0;C_0];

[t_2,y_2] = ode45(@(t,y) TumorCarOv(t,y,par),time_2,Init_2,[]);

A_num=0;
t_area=[t_1;t_2(2:end)];
y_area=[y_1(:,1:2);y_2(2:end,1:2)];
for s=1:(size(t_area,1)-1)
    h=(t_area(s+1)-t_area(s));
    b1=y_area(s,1)+y_area(s,2);
    b2=y_area(s+1,1)+y_area(s+1,2);
    A_num=A_num+(b1+b2)*h/2;
end
Area_comb(i+1,k+1)=A_num;
FinalT_comb(i+1,k+1)=y_area(end,1)+y_area(end,2);
end

else 
if adm_time_C>1
time_1=T_point(1:adm_time_C-1);
T_0=PBT030_ave(1);
C_0=0;
I_0=0;
V_0=0;
Init_1=[T_0;I_0;V_0;C_0];

[t_1,y_1] = ode45(@(t,y) TumorCarOv(t,y,par),time_1,Init_1,[]);

time_2=T_point(adm_time_C-1:end);
T_0=PBT030_ave(1);
C_0_cellNum=((T_0-0.1602)/0.0001946)/25;
C_0=C_0_cellNum*0.0001946+0.1602;
T_0=y_1(end,1);
I_0=0;
V_0=0.002;
Init_2=[T_0;I_0;V_0;C_0];

[t_2,y_2] = ode45(@(t,y) TumorCarOv(t,y,par),time_2,Init_2,[]);

A_num=0;
t_area=[t_1;t_2(2:end)];
y_area=[y_1(:,1:2);y_2(2:end,1:2)];
for s=1:(size(t_area,1)-1)
    h=(t_area(s+1)-t_area(s));
    b1=y_area(s,1)+y_area(s,2);
    b2=y_area(s+1,1)+y_area(s+1,2);
    A_num=A_num+(b1+b2)*h/2;
end
Area_comb(i+1,k+1)=A_num;
FinalT_comb(i+1,k+1)=y_area(end,1)+y_area(end,2);

else

time_1=T_point(1:end);
T_0=PBT030_ave(1);
C_0_cellNum=((T_0-0.1602)/0.0001946)/25;
C_0=C_0_cellNum*0.0001946+0.1602;
I_0=0;
V_0=0.002;
Init_1=[T_0;I_0;V_0;C_0];

[t_1,y_1] = ode45(@(t,y) TumorCarOv(t,y,par),time_1,Init_1,[]);

A_num=0;
t_area=t_1;
y_area=y_1(:,1:2);
for s=1:(size(t_area,1)-1)
    h=(t_area(s+1)-t_area(s));
    b1=y_area(s,1)+y_area(s,2);
    b2=y_area(s+1,1)+y_area(s+1,2);
    A_num=A_num+(b1+b2)*h/2;
end
Area_comb(i+1,k+1)=A_num;
FinalT_comb(i+1,k+1)=y_area(end,1)+y_area(end,2);
end 

end

end
end

% Save the results of the simulation to a MAT-file
save('Timing_CAR50_OV002.mat', 'X', 'Y', 'Area_comb', 'FinalT_comb')

function [dydt] = TumorCarOv(t, y, par)
    % Define the global variable b for CAR T-cell effector concentration.
    global b
    
    % Extract model parameters from the input 'par' array.
    alfa = par(1);      % Growth rate of tumor cells
    K = par(2);         % Carrying capacity (maximum tumor size)
    beta = par(3);      % Infection rate of tumor cells by virus
    theta1 = par(4);    % CAR T-cells killing of tumor 
    gamma = par(5);     % Tumor death rate 
    theta2 = par(6);    % CAR T-cells killing of infected tumor 
    b = par(7);         % Burst size
    theta3 = par(8);    % Infection rate of CAR T-cells 
    omega = par(9);     % Decay rate of the virus
    delta = par(10);    % Death rate of CAR T-cells
    theta4 = par(11);   % Expansion/exhaustion rate of CAR T cells
    
    % Define state variables:
    T = y(1); % Tumor population
    I = y(2); % Immune cells (CAR T-cells)
    V = y(3); % Virus concentration
    C = y(4); % CAR T-cell concentration
    
    % Differential equations:
    
    % Rate of change of tumor size (T):
    dT = alfa * T * (1 - (T + I) / K) - beta * T * V - theta1 * C * T;
    
    % Rate of change of immune (CAR T) cell population (I):
    dI = beta * T * V - gamma * I - theta2 * I * C;
    
    % Rate of change of virus concentration (V):
    dV = gamma * b * I - omega * V;
    
    % Rate of change of CAR T-cell concentration (C):
    dC = -delta * C + theta4 * C * (T + I) - theta3 * C * V;
    
    % Return the system of differential equations as a column vector:
    dydt = [dT; dI; dV; dC];
end