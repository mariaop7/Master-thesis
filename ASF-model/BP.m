function W = BP(xpos, kx, weights)
% Function for calculating the beampattern of an array
% IN:
%       xpos        Position of array elements
%       kx          Wavenumber samples
%       weights     Weights on the array elements
% OUT:
%       W           Beampattern of the array

% Ensuring all vectors are column vectors
xpos = xpos(:);
weights = weights(:); 
kx = kx(:);

m = size(kx,1); % Nr. of samples in the beampattern
W = zeros(1,m);

% Calculate beampattern
for i=1:m
    W(i) = sum(weights.*exp(1j*kx(i)*xpos));
end
