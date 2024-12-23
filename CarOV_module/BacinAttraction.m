%% Basin of attraction in timing plots

% Clear workspace and close any open figures
clear all
close all

% Load required data files
load('Timing_CAR25_OV002.mat')
load('Areas_CAR25_OV002.mat')

% Find the maximum and minimum values from the experimental data
Max = max([A_exp_Y1, A_exp_Y2, A_exp_Y3]);  % Maximum value across all experimental data
Min = min([A_exp_Y1, A_exp_Y2, A_exp_Y3]);  % Minimum value across all experimental data

% Initialize index counter for valid points
k = 1;

% Loop through the Area_comb matrix to filter values within the Min-Max range
for i = 1:size(Area_comb, 1)  % Iterate over rows
    for j = 1:size(Area_comb, 2)  % Iterate over columns
        % Check if the value in Area_comb is within the desired range
        if Area_comb(i, j) <= Max && Area_comb(i, j) >= Min
            x(k) = i;  % Store row index
            y(k) = j;  % Store column index
            k = k + 1;  % Increment counter for valid points
        end
    end
end

% Plot the valid (x, y) points as a scatter plot
figure
plot(x, y, '.', 'MarkerSize', 10)
hold on

% Compute and plot the convex boundary around the points
Bacin = boundary(x', y');  % Find the convex hull of the (x, y) points
plot(x(Bacin), y(Bacin), 'r-', 'LineWidth', 2)  % Plot the boundary in red


% Adjust plot appearance
axis equal  % Equal scaling on both axes
grid on  % Enable grid for better readability