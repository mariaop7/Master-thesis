function rel_power = rel_blob_power(img,N_rx, y0, fc,c,bw,lambda,range, u,d, save, draw, avg)
% Calculates the relative sum of power or average power for mainlobe, grating lobe 1, etc... 
% IN:
%       img (u,range)(complex array)    : image you want to find the relative blob power 
%       save                            : if save plots, save = 1
%       draw                            : if superpose blobs on image, draw = 1
%       avg                             : whether calculate average power over blobs, avg = 1, otherwise sum.
% OUT:
%       rel_power (array)               : relative power for grating lobe 1, 2, 3 etc.  

[U, R] = meshgrid(u, range);
nr_globes = round(d/lambda) - 1;
if d/lambda < 1.5
    nr_globes = 1;
end
L = N_rx*d-d;
flipimg = img.';
[Nu,Nr] = size(img);
globe_idx = [-nr_globes:-1, 1:nr_globes];
blobs = zeros(Nr,Nu, nr_globes*2);
angle_globew = 2*lambda/(N_rx*d);

if draw 
    % Plot fig
    f=figure('Position', [360 198 767.3333 406]);
    img_db = db(abs(img))-max(db(abs(img)),[], 'all');
    clims = [-50 0];
    imagesc(u, range, img_db.', clims);
    hold on
    
    L = N_rx*d-d;
    L_lmb = L/lambda;
    limit_lmb = ((N_rx-1)*d/lambda)^2;
    y_rel = y0/(L^2/lambda);
    rd = c/(2*bw);
    
    cb = colorbar(); 
    cb.Label.VerticalAlignment = "bottom";
    cb.LineWidth = 1;
    ylabel(cb,'Power [dB]','FontSize',11,'Rotation',270, 'FontName', 'Serif')
    title('Beamformed image')
    subtitle(sprintf('B/fc = $%.2f$, %.f elements, $d/\\lambda = %.1f$, $L = %.f\\lambda$, $L/\\Delta r = %.2f$, $L^2/\\lambda = %.f\\lambda$, $\\frac{y}{L^2/\\lambda} = %.1f$', bw/fc, N_rx,d/lambda, L_lmb, L/rd,  limit_lmb, y_rel),'Interpreter','latex')
    
    xlabel('$u$', 'Interpreter', 'latex')
    ylabel('Range [m]')
    set(gca,'LineWidth', 1, 'FontName', 'Serif', 'Fontsize', 11)
end

sets = 0;
for n=1:length(globe_idx)
    % Iterate over blobs from leftmost to rightmost (excluding mainlobe)

    if c/(2*bw) > abs(globe_idx(n))*(c/(fc+bw/2))/d*L/2
        % If range spread is less than range resolution, use range
        % resolution
        set1 = (R > y0-c/(2*bw));
        set2 = (R < y0+c/(2*bw));
    else
        set1 = (R > L/4*U + y0); 
        set2 = (R < -L/4*U + y0);
    end

    set3 = set1 & set2;
    set4 = -set1==0 & -set2==0;
    set5 = set3+set4; % Range spread 

    if angle_globew < abs(globe_idx(n))*((c/(fc-bw/2))/d - (c/(fc+bw/2))/d)
        % If the width due to broadband is larger than mainlobe width 
        set6 = U > globe_idx(n)*c/(fc+sign(globe_idx(n))*bw/2)/d;
        set7 = ~(U > globe_idx(n)*c/(fc-sign(globe_idx(n))*bw/2)/d);
    else
        set6 = U > globe_idx(n)*(c/fc)/d-(c/fc)/(N_rx*d);
        set7 = ~(U > globe_idx(n)*(c/fc)/d+(c/fc)/(N_rx*d));
    end

    set8 = set6 & set7;
    set9 = set5 & set8; % Blob indices

    sets = sets + set9;
    blobs(:,:, n) = set9; % Save indices

    if draw
        % Draw trapezoids on top
        if angle_globew < abs(globe_idx(n))*((c/(fc-bw/2))/d - (c/(fc+bw/2))/d)
            xleft = globe_idx(n)*(c/(fc+sign(globe_idx(n))*bw/2))/d;
            xright = globe_idx(n)*(c/(fc-sign(globe_idx(n))*bw/2)/d);
        else
            xleft = globe_idx(n)*lambda/d-lambda/(N_rx*d);
            xright = globe_idx(n)*lambda/d+lambda/(N_rx*d);
        end

        y1 = median(range)-L/4*xleft;
        y2 = median(range)-L/4*xright;
        y3 = median(range)+L/4*xright;
        y4 = median(range)+L/4*xleft; 
        trapes_x = [xleft, xright, xright, xleft, xleft];
        trapes_y = [y1, y2, y3, y4, y1];
        plot(trapes_x, trapes_y, 'r','LineWidth', 1)

        mainlobe_box = [-lambda/(N_rx*d) y0-c/(2*bw) 2*lambda/(N_rx*d) c/bw];
        rectangle('Position', mainlobe_box,'EdgeColor', 'r', 'LineWidth', 1)
    end
end

mainlobe_us = (abs(u)<lambda/(N_rx*d));
mainlobe_range = (range<median(range)+c/(2*bw) & range>median(range)-c/(2*bw));
x_u = find(mainlobe_us);
y_r = find(mainlobe_range);
main_blob = img(x_u(1):x_u(end), y_r(1):y_r(end));

if avg
    % Returns the relative average power
    avg_amp = zeros(1,nr_globes+1); % From largest grating lobe blob -> mainlobe
    for n=1:nr_globes
        blob = abs(flipimg(logical(blobs(:,:,n)))).^2;
        avg_amp(n) = mean(blob,'all'); % Mean over amplitude values
    end

    avg_amp(end) = mean(abs(main_blob).^2,'all');
%     rel_power = zeros(1,nr_globes);
% 
%     % Calculate relative power
%     for n=1:nr_globes
%         rel_power(n) = db(avg_amp(n)/avg_amp(end));
%     end
    rel_power = flip(avg_amp); % Mainlobe, 1st grating lobe, 2nd grating lobe...

else
    % Returns the sum power inside boxes
    sum_power = zeros(1,nr_globes+1); % From largest grating lobe blob -> mainlobe
    for n=1:nr_globes
        blob = abs(flipimg(logical(blobs(:,:,n)))).^2; % Power
        sum_power(n) = sum(blob,'all');
    end
    
    sum_power(end) = sum(abs(main_blob).^2, 'all');
    rel_power = flip(sum_power);  % main lobe, 1st grating lobe, 2nd grating lobe
end

if save
    % Whether save image
    figname = sprintf('./Code/figs/grating-lobe-analysis/4/imagefullwblobs_rb%.2f_y%.f_L%.1f_fc%.f_d%.1f.png', bw/fc, y0, N_rx*d-d, fc, d/lambda);
    exportgraphics(f, figname, 'Resolution', 300);
    figname = sprintf('./Code/figs/grating-lobe-analysis/4/imagefullwblobs_rb%.2f_y%.f_L%.1f_fc%.f_d%.1f.pdf', bw/fc, y0, N_rx*d-d, fc, d/lambda);
    exportgraphics(f, figname, 'Resolution', 300);
end