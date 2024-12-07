%-*-*-*-*-*-*-*-*-*-*-*-*-*-*
% Part 3: Least Squares
% Isela Gonzales
%-*-*-*-*-*-*-*-*-*-*-*-*-*-*

clc;

step = 1; % Step size for model
startDay = 0; % start of simulation
endDay = 30; % end pf simulation
time = startDay:step:endDay; % time vector corresponding to days
beta_t = 0.3; % True transimission rate
gamma_t = 0.1; % True recovery rate

%Initial conditions from part 1%
totalPopulation = 1000; % Total # of people
S0 = 990; % Initial # of susceptible people
I0_t = 10; % True inital # of infected people
R0 = 0; % Initial # of recovered people

% Function from part 1
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

% Initialize Variables
S = S0;
R = R0;
I = I0_t;
I_t = zeros(length(time));
I_t(1) = I0_t;

% Loop for non-linear model
for i = 2:length(time)
    [dS, dI, dR] = sir_derivatives(S, I, R, beta_t, gamma_t, totalPopulation);
    S = S + step * dS;
    I = I + step * dI;
    R = R + step * dR;
    I_t(i) = I;
end

I_lin = log(I_t); % Linearize I(t)

% Least squares application
x = [ones(length(time), 1), time'];
y = I_lin;
a = (x \ y); % coefficient value for regression
I0_linEst = a(1); % Estimate of ln(I0)
kEstimate = a(2); % k estimate
I0Estimate = exp(I0_linEst); % I0 value

% Estimation of beta and k values
k_t = (beta_t * S0 / totalPopulation) - gamma_t;
betaEstimate = (kEstimate + gamma_t) * totalPopulation / S0;

% Display results
disp('Data for 30 days: ')
disp('Estimated Beta: ');
disp(betaEstimate);
disp('True Beta: ');
disp(beta_t);
disp('Estimated k value: ');
disp(kEstimate);
disp('True k value: ');
disp(k_t);
disp('Estimated I0: ');
disp(I0Estimate);
disp('True I0: ');
disp(I0_t);

% Repeat previous calculations for 10 days :)
time_new = time(1:10); % using first 10 values of time matrix
I_linNew = I_lin(1:10);

% Least squares application for 10 days
x_new = [ones(length(time_new), 1), time_new']; %new x values
y_new = I_linNew'; %new y values
a_new = (x_new \ y_new);
I0_linEst_new = a_new(1);
kEstimate_new = a_new(2);

I0Estimate_new = exp(I0_linEst_new); % new I0 estimate
betaEstimate_new = (kEstimate_new + gamma_t) * totalPopulation / S0;

% Display results
disp('Data for 10 days: ')
disp('Estimated Beta: ');
disp(betaEstimate_new);
disp('True Beta: ');
disp(beta_t);
disp('Estimated k value: ');
disp(kEstimate_new);
disp('True k value: ');
disp(k_t);
disp('Estimated I0: ');
disp(I0Estimate_new);
disp('True I0: ');
disp(I0_t);

% Discussion:
% Based on the results, the data from merely 10 days of running the model
% deemed to be more accurate in comparison to 30 days of data. This could
% be due to the values needed being at the beginning stages of the model
% and the linearization of this model makes the more accurate numbers in
% the beginning.