function omega_seg = abel_smith_divide(alpha,N_IIR)
% function calculates the segmentation for the Abel-Smith algorithm
% alpha: channel characteristic alpha=lambda_0^2*B^2*D*L/(4*pi*c)
% N_IIR: number of first order IIR all-pass sections
% omega_seg = [omega_0,omega_1,...,omega_N_IIR]: vector of frequencies omega_i separating the 2*pi-area frequency bands
omega_seg(1,1) = -pi;
b = -N_IIR;
a = alpha;
for i = 2:N_IIR+1
    w0 = omega_seg(1,i-1);
    c = 2*pi + w0*N_IIR-alpha*w0^2;
    inside = (b^2)-(4*a*c);
    omega_seg(1,i) = (-b-sqrt(inside))/(2*a);
end

end

