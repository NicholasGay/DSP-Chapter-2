function bit_error_rate(h_tot,index,rho,theta,phi)

%% Simulation parameters
%% QPSK, modulation Format
M = 4; 
hMod = modem.pskmod('M', M, 'PhaseOffset', pi/4, 'SymbolOrder', 'Gray','InputType', 'Bit');

% QPSK, demodulation format
hDemod = modem.pskdemod('M', M, 'PhaseOffset', pi/4,...
    'SymbolOrder', 'Gray', 'OutputType', 'Bit');

%% Simulation
SNR=-4:2:14;                % in dB
len_SNR=length(SNR);

% Original message
numSymb = 4*25000;  % Number of symbols
msg_orig = (sign(randn(numSymb,1))+1)/2; 

% QPSK modulation
 out_qpsk = modulate(hMod,msg_orig);
        
% Signal through the channel
out_channel=conv(h_tot,out_qpsk);        

% AWG Noise
channel_noise = (randn(length(out_channel),1)+1i*randn(length(out_channel),1))/sqrt(2); 

for j=1:len_SNR  
                
        % Channel + Noise    
        sigma_n=1/sqrt(10^(SNR(j)/10)); % variance of noise
        out_channel_plus_noise=out_channel+channel_noise*sigma_n;
        
        % Convolving output channel plus noise with the all-pass Equalizer
        % Task 3: Students just need to use their code from Task 3 and
        % change the input
        output_equalizer=conv_anyinput_allpass_equalizer(rho,theta,phi,out_channel_plus_noise);
                
        % Parsing out the desired signal which was send before the CD channel
        output_equalizer=output_equalizer(index:index+length(out_qpsk)-1);
                
        % QPSK demodulation
        msg_demod = demodulate(hDemod, output_equalizer);
                
        % Bit Error Rate
        [errorBit ratioBit(j)] = biterr(msg_orig, msg_demod);        

end

figure
semilogy(SNR,ratioBit)
hold on
semilogy(SNR,0.5*erfc(sqrt(10.^(0.1*SNR)/2)),'--r')
grid
title('BER vs SNR')
xlabel('SNR (dB)')
ylabel('BER')
legend('IIR Equalizer','AWGN channel')

end
