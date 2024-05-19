function [area, a, b, h] = trapezoid_area(N_rx, bw, fc, c, u, d)
% Calculates the area of a trapezoid of grating lobe at u.

lambda = c/fc;
lambda_max = c./(fc-bw/2);
lambda_min = c./(fc+bw/2);

% p = d*u/lambda;
L = (N_rx-1).*d;
a = abs(lambda_min*L/2/lambda.*u); 
b = abs(lambda_max*L/2/lambda.*u); 
h = abs(u./lambda*(lambda_max-lambda_min)); 

a(a < c/bw) = c/bw; % If the smallest parallel side is smaller than range resolution
b(b < c/bw) = c/bw; % If the largest parallel side is smaller than range resolution
% h(h < lambda./(N_rx.*d)*2) = lambda./(N_rx(h < lambda./(N_rx.*d)*2).*d(h < lambda./(N_rx.*d)*2))*2;
ml_width = (h < lambda./(N_rx.*d)*2);
h(ml_width) = lambda./(N_rx.*d)*2; % If the width due to broadband is smaller than the mainlobe width


area = abs((a+b).*h/2);

% area = (N_rx-1).*(lambda_min+lambda_max).*u/lambda.*(lambda_max-lambda_min)/2;
% OR
% area = (lambda_min+lambda_max)*L/2*u.^2/lambda^2*(lambda_max-lambda_min)/2;