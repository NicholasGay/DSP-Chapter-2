function y = MSE_TP(x,alpha,beta,N_IIR)
% function computes MSE in transfer function MSE_trans. phase
% x = [rho;theta;phi]
% rho:   column vector of N_IIR radii rho_i of poles of CD compensation filter
% theta: column vector of N_IIR angles theta_i of poles of CD compensation filter
% phi:   phase correction term phi_0 of CD compensation filter
% alpha: parameter alpha=lambda_0^2*B^2*D*L/(4*pi*c) characterizing CD channel
% beta:  integer factor beta=ceil(2*alpha*pi)
% N_IIR: number of first order IIR all-pass sections
% y:     MSE in transfer function MSE_trans. phase

w = linspace(-pi,pi,(2^14));
for i = 1:length(w)
    b0 = ones(N_IIR,1);
    z_inv(1:N_IIR,1) = exp(-1i*w(1,i));
    
    alpha0 = x(1:N_IIR,:).*exp(1i.*x(N_IIR+1:end-1,:));
    num = -conj(alpha0)+z_inv;
    dem = b0-alpha0.*z_inv;
    
    Tau_ideal = exp(-1i.*(x(end,1)+beta*w(1,i)));
    
    H = exp(-1i*alpha*w(1,i)^2);
    GH = prod(num./dem)*H;
    
    Y(1,i)=abs(GH-Tau_ideal)^2;
end
y = trapz(Y);
end

