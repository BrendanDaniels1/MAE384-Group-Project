%==========================================
% Variable Declaration
%==========================================

clear;clc;

% Array of disease names for reference during plotting
diseases = {'Influenza', 'COVID-19', 'Measles'};

% Parameters for each disease: [Beta, Gamma]
% Beta: transmission rate, Gamma: recovery rate
parameters = [0.3, 0.1;  % Influenza
              1.0, 0.1;  % COVID-19
              2.0, 0.2]; % Measles

% Initial conditions for the population
totalPopulation = 1000; % Total number of people in the population
S0 = 990; % Initial number of susceptible individuals
I0 = 10;  % Initial number of infected individuals
R0 = 0;   % Initial number of recovered individuals

% Simulation time parameters
startDay = 0;    % Start of the simulation (day 0)
endDay = 100;    % End of the simulation (day 100)
step = 1;        % Time step size (1 day)

% Time vector representing each simulation step
time = startDay:step:endDay;

%==========================================
% Main Loop for Each Disease
%==========================================

for i = 1:3 % Loop through the three diseases
    % Extract the parameters Beta (transmission rate) and Gamma (recovery rate)
    Beta = parameters(i, 1);
    Gamma = parameters(i, 2);

    % Initialize the state variables for SIR
    Sx = S0; % Susceptible population
    Ix = I0; % Infected population
    Rx = R0; % Recovered population

    % Create arrays to store the S, I, and R values over time
    S_vals = zeros(1, length(time)); % Susceptible
    I_vals = zeros(1, length(time)); % Infected
    R_vals = zeros(1, length(time)); % Recovered

    % Assign initial conditions to the first element of the arrays
    S_vals(1) = Sx;
    I_vals(1) = Ix;
    R_vals(1) = Rx;

    %--------------------------------------------
    % Numerical Simulation Using RK4
    %--------------------------------------------
    for t = 2:length(time) % Loop over time steps (from the second step onwards)
        % Calculate RK4 updates for S, I, and R
        [k1_S, k1_I, k1_R] = sir_derivatives(Sx, Ix, Rx, Beta, Gamma, totalPopulation);
        [k2_S, k2_I, k2_R] = sir_derivatives(Sx + 0.5 * step * k1_S, ...
                                             Ix + 0.5 * step * k1_I, ...
                                             Rx + 0.5 * step * k1_R, ...
                                             Beta, Gamma, totalPopulation);
        [k3_S, k3_I, k3_R] = sir_derivatives(Sx + 0.5 * step * k2_S, ...
                                             Ix + 0.5 * step * k2_I, ...
                                             Rx + 0.5 * step * k2_R, ...
                                             Beta, Gamma, totalPopulation);
        [k4_S, k4_I, k4_R] = sir_derivatives(Sx + step * k3_S, ...
                                             Ix + step * k3_I, ...
                                             Rx + step * k3_R, ...
                                             Beta, Gamma, totalPopulation);

        % Update S, I, and R using the RK4 weighted average formula
        Sx = Sx + (step / 6) * (k1_S + 2 * k2_S + 2 * k3_S + k4_S);
        Ix = Ix + (step / 6) * (k1_I + 2 * k2_I + 2 * k3_I + k4_I);
        Rx = Rx + (step / 6) * (k1_R + 2 * k2_R + 2 * k3_R + k4_R);

        % Store the updated values into the arrays
        S_vals(t) = Sx;
        I_vals(t) = Ix;
        R_vals(t) = Rx;
    end

    %--------------------------------------------
    % Plot Results for Current Disease
    %--------------------------------------------
    figure; % Create a new figure for each disease
    plot(time, S_vals, 'b', 'LineWidth', 2); % Plot Susceptible population in blue
    hold on; % Allow multiple plots on the same figure
    plot(time, I_vals, 'r', 'LineWidth', 2); % Plot Infected population in red
    plot(time, R_vals, 'g', 'LineWidth', 2); % Plot Recovered population in green
    title(['Disease Progression: ', diseases{i}]); % Set title with disease name
    xlabel('Time (days)'); % Label x-axis
    ylabel('Population'); % Label y-axis
    legend('Susceptible', 'Infected', 'Recovered'); % Add legend
    grid on; % Enable grid for better visualization
end

%==========================================
% Function to Compute Derivatives
%==========================================

function [dS, dI, dR] = sir_derivatives(S, I, R, Beta, Gamma, N)
    % This function calculates the derivatives dS/dt, dI/dt, and dR/dt
    % Inputs:
    %   S - Susceptible population
    %   I - Infected population
    %   R - Recovered population
    %   Beta - Transmission rate
    %   Gamma - Recovery rate
    %   N - Total population (constant)

    dS = -Beta * S * I / N;      % Rate of change of susceptible population
    dI = Beta * S * I / N - Gamma * I; % Rate of change of infected population
    dR = Gamma * I;              % Rate of change of recovered population
end
