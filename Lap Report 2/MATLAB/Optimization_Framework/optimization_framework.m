function [rho,theta,phi] = optimization_framework(alpha)
%% Implement the optimization framework here
%  alpha: parameter alpha=lambda_0^2*B^2*D*L/(4*pi*c) characterizing CD channel
%  rho:   column vector of N_IIR radii rho_i of poles of CD compensation filter
%  theta: column vector of N_IIR angles theta_i of poles of CD compensation filter
%  phi:   phase correction term phi_0 of CD compensation filter


% Calculate N_IIR, beta
N_IIR = ceil(2*pi*alpha);
beta = N_IIR;

%% 1) Abel-Smith algorithm for computing intial rho_i and theta_i

% Uncomment the following function
[rho,theta] = abel_smith(alpha,N_IIR);



%% 2) Non-linear optimization of MSE in group delay for refining rho_i and theta_i
options = optimoptions('fmincon','Algorithm','active-set');
x_init = [rho;theta];
lb = zeros(size(x_init));
lb(1:N_IIR,1,1) = 0;
lb(N_IIR+1:end,1) = -pi;

ub = ones(size(x_init));
ub(1:N_IIR,1) = 0.98;
ub(N_IIR+1:end,1) = pi;
x_opt = fmincon(@(x)(MSE_GD(x,alpha,beta,N_IIR)),x_init,[],[],[],[], lb,ub,[],options);


%% 3) Non-linear optimization of MSE in transfer function for computing initial phi_0
%phi = 0;
%y = MSE_TP([x_opt;phi],alpha,beta,N_IIR);

phi_guess = pi/2;

phi_init = fmincon(@(x)(MSE_TP([x_opt;x],alpha,beta,N_IIR)),phi_guess,[],[],[],[],-pi,pi,[],options);

%% 4) Non-linear optimization of MSE in transfer function for computing optimal rho_i, theta_i and phi_0
x3_init= [x_opt;phi_init];
lb = zeros(size(x3_init));
lb(1:N_IIR,1) = 0;
lb(N_IIR+1:end,1) = -pi;

ub = ones(size(x3_init));
ub(1:N_IIR,1) = 0.98;
ub(N_IIR+1:end,1) = pi;

x_opt_final = fmincon(@(x)(MSE_TP(x,alpha,beta,N_IIR)),x3_init,[],[],[],[],lb,ub,[],options);

rho = x_opt_final(1:N_IIR,1);
theta = x_opt_final(N_IIR+1:end-1,1);
phi = x_opt_final(end,1);

end