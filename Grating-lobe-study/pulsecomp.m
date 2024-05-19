function pc_data = pulsecomp(x, s)
% x : data
% s : transmit signal

[pc_data, lags] = xcorr(x,s);

% Pick positive lags 
idx = lags>=0;
pc_data = pc_data(idx);

end

