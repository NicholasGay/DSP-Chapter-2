function main_program()

clc
addpath('Optimization_Framework','CD_Compensation')

%% Operating wavelength and related constants
lambda=1550*1e-9;                      % Operating wavelength in meters
c=3e8;                                 % Speed of light in m/s
D=16*1e-12/(1e-9*1e3);                 % Dispersion in ps/(nm*km)

%% Configuration Parameters For the CD Channel
L=23*1e3;                              % Length of fiber in meters
Bw=56e9;                               % Bandwidth of the signal in Hz

%% Sampling frequency
B=Bw;

%% Calculating alpha
alpha=lambda^2*B^2*D*L/(4*pi*c);        % Multiplying the above parameters

%% Task 1: (For students)
%  Implement the optimization framework for computing the optimal
%  parameters of the CD compensation filter
%% Computing optimal parameters of CD compensation filter
[rho,theta,phi] = optimization_framework(alpha);

%% Number of frequency points needed to generate channel impulse response
freq_pts=1024;

%% Task 2: (For students)
%  Write a function which generates the impulse response of the CD channel
%% Generating channel impulse response
h_CD=impulse_response_channel(alpha,freq_pts); 
figure
stem(0:(freq_pts-1),abs(h_CD))
grid on
title('Impulse response of CD channel')
xlabel('n')
ylabel('|h_{CD}[n]|')

%% Task 3: (For students)
%  Write a function which convolves any input with all-pass equalizer
%% In this case, input is the impulse response of CD channel
GH=conv_anyinput_allpass_equalizer(rho,theta,phi,h_CD);
[val,index]=max(abs(GH));

%% Plot GH, ideally it should be an impulse
figure
stem(0:(freq_pts-1),abs(GH))
grid on
title('Output of IIR equalizer')
xlabel('n')
ylabel('|gh[n]|')

%% BER
bit_error_rate(h_CD,index,rho,theta,phi)

end
