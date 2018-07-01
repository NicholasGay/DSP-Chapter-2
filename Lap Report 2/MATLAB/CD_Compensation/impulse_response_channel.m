function h_tot = impulse_response_channel(alpha,freq_pts)
%% Generate the impulse response in time domain
% alpha: channel characteristic alpha=lambda_0^2*B^2*D*L/(4*pi*c)
% freq_pts: number of frequency points
% h_tot: column vector with time domain values

w = linspace(-pi,pi,freq_pts);
%% Sampling Frequency Response
H_CD = exp(-1i*alpha*(w.^2));

%% IFFT and normalization

h_tot = ifft(H_CD).';
% return column vector


end