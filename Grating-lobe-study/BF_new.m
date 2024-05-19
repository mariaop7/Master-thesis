function [image_full] = BF_new(mfdata, x_grid, y_grid, x_pos, x_tx, fs, N_rx, c, M, range, theta)

n_r = length(range);
n_theta = length(theta);

dt = 1/fs;
taxe = 0:dt:M*dt-dt;

image_full = complex(zeros(n_theta,n_r));
for n=1:N_rx % Number of receivers
    tmp_img = complex(zeros(n_theta,n_r));
    for th=1:n_theta
        r_t = sqrt( (x_tx-x_grid(th, :)).^2 + (y_grid(th, :)).^2 );
        r_r = sqrt( (x_pos(n)-x_grid(th, :)).^2 + (y_grid(th, :)).^2 );
        r = r_t + r_r; t_delay = r / c;
        tmp_img(th,:) = interp1(taxe, mfdata(:,n), t_delay); % Linear interpolation
    end
    image_full = image_full + tmp_img;
end % n rx
