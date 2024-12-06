%Define Beta Function for transmission rate
T_r = .3.*(1+5.*sin(2.*pi.*(365/365).*time));
Beta(i) = .3;
Gamma = .1 ;
h = .1; %step size
t0 = 0;
tf = 30 ;
time = t0:h:tf;
%Population initial conditions
Total = 1000;
S0 = 990;
I0 = 10;
R = 0;
%Arrays
S_1 = zeros(1,length(time));
I_1= zeros(1,length(time));
R_1= zeros(1,length(time));
S_1(1)= S0;
I_1(1)= I0;
R_1(1)= R0;

for i = 2:length(time)
    [Sk1,Ik1, Rk1] = dSIRdt(Sx, Ix, Rx,Beta(i), Gamma, Total);
    [Sk2,Ik2, Rk2] = dSIRdt(Sx + (.5*Sk1*h),Ix + (.5*Ik1*h),Rx+(.5*Rk1*h),Beta(i), Gamma, Total);
    [Sk3,Ik3, Rk3] = dSIRdt(Sx + (.5*Sk2*h),Ix + (.5*Ik2*h),Rx+(.5*Rk2*h),Beta(i), Gamma, Total);
    [Sk4,Ik4, Rk4] = dSIRdt(Sx + (Sk3*h),Ix + (Ik3*h),Rx+(Rk3*h),Beta(i), Gamma, Total);
% Runge-Katta Formula
    Sx = Sx + ((h/6)*(Sk1+2*Sk2+2*Sk3+Sk4));
    Ix = Ix + ((h/6)*(Ik1+2*Ik2+2*Ik3+Ik4));
    Rx = Rx + ((h/6)*(Rk1+2*Rk2+2*Rk3+Rk4));
%Storing Values
    S_1(i) = Sx;
    I_1(i) = Ix;
    R_1(i) = Rx;
end

% Creating plot for Beta Function

figure;
plot (time, S_1, 'r', time, I_1,'g',R_1,'b')
legend('Susceptible','Infected','Recovered')
xlabel('Time')
ylabel('Respective Populations')
title('Variable Spread Rate Omega = 2pi*365/365')
grid on
%Discussion
%Periodic fluctuation of function visible(high frequency)

% Define frequency range
t = 30;
N=length(time);
f = (1/t)*(0:(N/2));
% Calculate spectrum
spectrum=fft(I);
P2 = abs(spectrum/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
plot(f,P1);
title("Spectrum of Infected Cases, Omega = 2pi*365/365");
xlabel("frequency Coefficient")
ylabel("Absolute Value fft(I(t))")

%Discussion

clearvars;

%Part Two, omega = 2pi*100/365

%Define Beta Function for transmission rate
T_r = .3.*(1+5.*sin(2.*pi.*(100/365).*time));
Beta(i) = .3;
Gamma = .1 ;
h = .1; %step size
t0 = 0;
tf = 30 ;
time = t0:h:tf;
%Population initial conditions
Total = 1000;
S0 = 990;
I0 = 10;
R = 0;
%Arrays
S_1 = zeros(1,length(time));
I_1= zeros(1,length(time));
R_1= zeros(1,length(time));
S_1(1)= S0;
I_1(1)= I0;
R_1(1)= R0;

for i = 2:length(time)
    [Sk1,Ik1, Rk1] = dSIRdt(Sx, Ix, Rx,Beta(i), Gamma, Total);
    [Sk2,Ik2, Rk2] = dSIRdt(Sx + (.5*Sk1*h),Ix + (.5*Ik1*h),Rx+(.5*Rk1*h),Beta(i), Gamma, Total);
    [Sk3,Ik3, Rk3] = dSIRdt(Sx + (.5*Sk2*h),Ix + (.5*Ik2*h),Rx+(.5*Rk2*h),Beta(i), Gamma, Total);
    [Sk4,Ik4, Rk4] = dSIRdt(Sx + (Sk3*h),Ix + (Ik3*h),Rx+(Rk3*h),Beta(i), Gamma, Total);
% Runge-Katta Formula
    Sx = Sx + ((h/6)*(Sk1+2*Sk2+2*Sk3+Sk4));
    Ix = Ix + ((h/6)*(Ik1+2*Ik2+2*Ik3+Ik4));
    Rx = Rx + ((h/6)*(Rk1+2*Rk2+2*Rk3+Rk4));
%Storing Values
    S_1(i) = Sx;
    I_1(i) = Ix;
    R_1(i) = Rx;
end

% Creating plot for Beta Function

figure;
plot (time, S_1, 'r', time, I_1,'g',R_1,'b')
legend('Susceptible','Infected','Recovered')
xlabel('Time')
ylabel('Respective Populations')
title('Variable Spread Rate Omega = 2pi*365/365')
grid on
%Discussion
%Periodic fluctuation of function visible(high frequency)

% Define frequency range
t = 30;
N=length(time);
f = (1/t)*(0:(N/2));
% Calculate spectrum
spectrum=fft(I);
P2 = abs(spectrum/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
plot(f,P1);
title("Spectrum of Infected Cases, Omega = 2pi*100/365");
xlabel("frequency Coefficient")
ylabel("Absolute Value fft(I(t))")