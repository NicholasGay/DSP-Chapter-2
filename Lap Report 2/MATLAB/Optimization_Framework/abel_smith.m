function [rho,theta] = abel_smith(alpha,N_IIR)
% Abel-Smith algorithm
% alpha: channel characteristic alpha=lambda_0^2*B^2*D*L/(4*pi*c)
% N_IIR: number of first order IIR all-pass sections
% rho: column vector with N_IIR filter radii rho_i
% theta: column vector with N_IIR filter angles theta_i
% xi: facultative correction factor between 0.75 and 0.85

% Step 1: In order to find theta and rho, we need to calculate 
%         omega_seg = [omega_0,omega_1,...,omega_N_IIR] (see equations (21)
%         and (22))

omega_seg = abel_smith_divide(alpha,N_IIR); % Uncomment this function

% Step 2: Calculate theta and rho
sigma = 0.8;
for i = 1:N_IIR
    theta(i,1) = (omega_seg(1,i+1)+omega_seg(1,i))/2;
    delta = (omega_seg(1,i+1)-omega_seg(1,i))/2;
    niu = (1-sigma*cos(delta))/(1-sigma);
    rho(i,1) = niu - sqrt(niu^2-1);
end

end

