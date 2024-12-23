%% Bifurcation Diagram for CAR Module

% This script generates a bifurcation diagram by plotting solutions of the
% CAR T-cell system for varying values of parameter 'B'.

clear all;  % Clear all variables
close all;  % Close any open figures

% Create a new figure window
figure;

% Set parameter 'a' for the model
A = 2;

% Define various values for parameter B and corresponding X values

% Bifurcation points and branches
% First vertical line (b1) from B=0 to B=a+10 at X=0
B1 = [0, A + 10];  
B10 = [0, 0];  % X = 0, horizontal line

% Second branch (b2) from B=a to B=a+10, with corresponding X values (b20)
B2 = linspace(A, A + 10, 50);  % B values from a to a+10
B20 = A ./ B2;  % X values are inversely proportional to B

% Third branch (b3) with a constant X value of 1, and B from 0 to a
B3 = [0, A];
B30 = [1, 1];  % X = 1, constant horizontal line

% Fourth branch (b4) with a constant X value of 1, and B from a to a+10
B4 = [A, A + 10];
B40 = [1, 1];  % X = 1, constant horizontal line

% Plot the bifurcation diagram
% Plot each of the branches with different line styles
plot(B1, B10, 'k--', ...        % Dashed black line for B from 0 to a+10, X=0
     B2, B20, 'k', ...          % Solid black line for second branch
     B3, B30, 'k', ...          % Solid black line for third branch
     B4, B40, 'k--', ...        % Dashed black line for fourth branch
     A, 1, 'ro', ...            % Red circle marker for the critical point at (a, 1)
     'LineWidth', 1, 'MarkerSize', 10); % Line width and marker size

% Set axis limits for the plot
axis([0 10 0 2]);

% Define custom x-tick labels
xticks([0 2 4 6 8 10]);
xticklabels({0, 'A', 4, 6, 8, 10});

% Define custom y-tick labels
yticks([0 .5 1 1.5 2]);

% Set axis labels with LaTeX formatting for better presentation
xlabel('$B$', 'Interpreter', 'latex', 'FontSize', 15);  % Label for B-axis
ylabel('$X$', 'Interpreter', 'latex', 'FontSize', 15);  % Label for X-axis

% Title for the plot (optional)
% title('Bifurcation Diagram of CAR Module', 'Interpreter', 'latex', 'FontSize', 16);