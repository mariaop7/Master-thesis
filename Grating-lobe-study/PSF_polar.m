function image_full = PSF_polar(options)
% Applies pulse-echo beamfoming to reconstruct point scatterer 
% OUT:
%       image_full (struct)

arguments 
    options.N_rx     (1,1) = 8          % Number of elements
    options.bw       (1,1) = 1e2        % Bandwidth
    options.t_p      (1,1) = 10e-2      % Pulse length      
    options.fc       (1,1) = 10e3       % Center frequency
    options.c        (1,1) = 1500       % Speed of sound
    options.fs       (1,1) = 200e3*6    % Sampling frequency
    options.d        (1,1) = 2*0.0150   % Interelement distance
    options.y0       (1,1) = 3          % Point y-position
    options.x0       (1,1) = 0          % Point x-position
    options.x_tx     (1,1) = 0          % Tx x-position
    options.Nu       (1,1) = 1000       % u samples
end

lambda = options.c/options.fc;
d = options.d; N_rx = options.N_rx;
x_pos = (1:N_rx)*d-(N_rx+1)*d/2;
bw = options.bw;
t_p = options.t_p; fc = options.fc;
fs = options.fs; y0 = options.y0;
x_tx = options.x_tx; x0 = options.x0;
c = options.c; 
L = abs(x_pos(1)-x_pos(end));
Nu = options.Nu;
drange = c/(2*bw); 

if y0 < 10
    range = 0.1:drange/50:(y0+10);
else
    range = (y0-10):drange/10:(y0+10);
end
u = linspace(-1,1,Nu); % theta from -pi/2 to pi/2, equispaced
theta = asin(u)+pi/2; % theta from 0 to pi, not equispaced

% convert theta and range to x and y grid
[x_grid, y_grid] = pol2cart(theta.',range);

dt = 1/fs;
a = bw/t_p;
t = 0:dt:t_p;
s = exp(1j*2*pi*((fc-bw/2)*t + a/2*t.^2)); % Transmit signal, LFM pulse

max_r = L + (y0+10)*2;
max_tsample = round(max_r/c*fs);
M = 2*max_tsample+length(t); % Time samples

taxe = 0:dt:M*dt-dt; % LFM pulse time axis

if ((max_tsample+length(t))>M)
    error('Error. Not enough samples M.')
end

% Generate data
rawdata = complex(zeros(M,N_rx)); % Output rawdata [time, Rx]

% s_zeropad = [s zeros(1,M-length(s))];
r_t = sqrt( (x_tx-x0)^2 + (y0)^2 );
for n=1:N_rx 
    r_r = sqrt( (x_pos(n)-x0)^2 + (y0)^2 );
    r = r_t + r_r; t_delay = r / c; 
    s_interpol = exp(1j*2*pi*((fc-bw/2)*(taxe-t_delay) + a/2*(taxe-t_delay).^2)); % get the right phase

    t_samples = round(t_delay*fs); % Nearest neighbor interpolator
    s_interpol(1:t_samples) = 0;
    s_interpol(t_samples+length(s)+1:end) = 0; % Delays s_interpol to sample t_sample

    rawdata(:,n) = s_interpol;
end

% Pulse compression
mfdata = complex(zeros(M,N_rx)); % [time, Rx]

for nr=1:N_rx
    mfdata(:,nr) = pulsecomp(rawdata(:,nr), s);
end

% Beamforming
image_full.data = BF_new(mfdata, x_grid, y_grid, x_pos, x_tx, fs, N_rx, c, M, range, theta); % Beamformed data

image_full.fnum = y0/L; % F-number
image_full.lim = L^2/lambda; % Far-field near-field limit
image_full.x_grid = x_grid;
image_full.y_grid = y_grid;
image_full.u = u;
image_full.range = range;
image_full.theta = theta;
image_full.bw = bw;
image_full.fc = fc;
image_full.y0 = y0;
image_full.L = L;