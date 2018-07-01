
B = 56*(10^9);
landa = 1550*(10^-9);
D = 16*(10^(-12+9-3));
L = 23*(10^3);
c = 299792458;


for k = 0:255
    w(1,k+1) = (2*pi*k)/256;
end


alpha = ((landa^2)*(B^2)*D*L)/(4*pi*c);
Hcd = exp(-1i*alpha*(w.^2));
r = real(ifft(Hcd));
im = imag(ifft(Hcd));
figure,

ax1 = subplot(2,1,1);
stem(ax1,r)
xlim([0 255])
title(ax1,'Real')
xlabel(ax1,'Samples[n]')

ax2 = subplot(2,1,2);
stem(ax2,im)
xlim([0 255])
title(ax2,'Imag')
xlabel(ax2,'Samples[n]')

arg = ((landa^2)/(2*c))*D*L*(B^2);
N = floor(arg)*2+1;