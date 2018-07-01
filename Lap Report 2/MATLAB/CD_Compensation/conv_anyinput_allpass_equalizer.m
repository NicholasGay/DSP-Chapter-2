function output = conv_anyinput_allpass_equalizer(rho,theta,phi,input)
%% Convolution of any filter specified by rho/theta/phi with input
% rho: column vector with N_IIR filter radii rho_i
% theta: column vector with N_IIR filter angles theta_i
% phi: phase correction term phi_0
% input: signal to be filtered
% output: filtered signal
% Hint: Use function filter()
alpha_i = rho.*exp(1i.*theta);
N = size(rho,1);
sos = zeros(2*N,6);
sos(1:N,1) = -conj(alpha_i);
sos(1:N,2) = ones(N,1);
sos(1:N,4) = ones(N,1);
sos(1:N,5) = -alpha_i;
sos(N+1:end-1,2) = 1;
sos(N+1:end,4) = ones(N,1);
sos(end,2) = exp(-1i*phi);
output = sosfilt(sos,input);

end