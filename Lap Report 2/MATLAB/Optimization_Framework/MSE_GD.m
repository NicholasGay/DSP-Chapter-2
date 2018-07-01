function y = MSE_GD(x,alpha,beta,N_IIR)
% function computes MSE in group delay MSE_GD
% x = [rho;theta]
% rho:   column vector of N_IIR radii rho_i of poles of CD compensation filter
% theta: column vector of N_IIR angles theta_i of poles of CD compensation filter
% alpha: parameter alpha=lambda_0^2*B^2*D*L/(4*pi*c) characterizing CD channel
% beta:  integer factor beta=ceil(2*alpha*pi)
% N_IIR: number of first order IIR all-pass sections
% y:     MSE in group delay MSE_GD

w = linspace(-pi,pi,(2^14));

for i = 1:length(w)
    num = 1-x(1:5,:).^2;
    dem = 1+ x(1:N_IIR,:).^2-2*x(1:N_IIR,:).*cos(w(1,i)-x(N_IIR+1:end,:));
    summation = sum(num./dem);
    tau_desired = -2*alpha*w(1,i)+beta;
    Y(1,i) = (abs(tau_desired - summation))^2;
end
y = trapz(Y);
end

