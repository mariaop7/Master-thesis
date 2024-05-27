function [total_BP, angles, u] = ASF_PW(options)
% Calculates the broadband aperture smoothing function (ASF)
% OUT:
%       total_BP        Total broadband ASF
%       angles          Corresponding angles in degrees
%       u               Corresponding u-values 
arguments 
    options.N           (1,1) = 28 % Receive elements
    options.M           (1,1) = 38 % Transmits
    options.bw          (1,1) = 3e4 % Bandwidth
    options.fc          (1,1) = 1e5 % Center frequency
    options.c           (1,1) = 1500 % Speed of sound
    options.d_rx        (1,1) = 2.5*(1500/1e5) % Receive element width/distance
    options.d_tx        (1,1) = 1.5*2.5*(1500/1e5) % Transmitter element width
    options.skew_angle  (1,1) = 0.41; % Yaw error
    options.L_tx        (1,1) = 0.525; % Distance between pings
    options.NFFT        (1,1) = 2048*6; % Number of samples in the ASF
    options.u_max       (1,1) = 1; % Maximum absolute u value in u-space
    options.w_Tx        (1,:) = hamming(1,38)/38; % Transmit weighitng
end

c = options.c; fc = options.fc;
bw = options.bw; 
N = options.N; M = options.M;
skew_angle = options.skew_angle;
L_tx = options.L_tx;
d_Rx = options.d_rx; d_Tx = options.d_tx;
u_max = options.u_max; NFFT = options.NFFT;
w_Tx = options.w_Tx;

f0 = fc - bw/2;
f1 = fc + bw/2;

freq = linspace(f0,f1,301);

% LFM spectrum
lfmspec = ones(1,length(freq)).*(tukeywin(length(freq), 0.1).'); % Tukey tapered
lfmspec = lfmspec.*(hamming(length(freq)).'); % Hamming tapered

total_BP = zeros(1,NFFT);
u = linspace(-u_max,u_max,NFFT);
angles = rad2deg(asin(u));

for f=1:length(freq) % Loop over frequencies 
   
    [W_2w, ~] = grating_lobe_BP(N, M, skew_angle, L_tx, d_Rx, d_Tx, c, freq(f),u, 0, w_Tx); % Single-frequency ASF

    pulse_factor = lfmspec(f); % Frequeny component of LFM pulse

    W_pulse = W_2w.*pulse_factor; % Multiply ASF by the frequency component

    total_BP = total_BP + W_pulse;
end