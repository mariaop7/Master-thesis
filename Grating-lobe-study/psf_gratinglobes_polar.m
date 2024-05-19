% Example use
bw = 3e3; fc = 1e4; 
c = 1500; t_p = 10e-2; N_rx = 8;
Nu = 500; fs = 200e3;
lambda = c/fc; d = lambda*2; 

% u = 0.2058; % U-position of first grating lobe
% d = lambda/u;
% N_rx = round(126*lambda/d+1);
L = N_rx*d-d;
% fnum = 2.5;
% y0 = fnum*L;
y0 = 50; 

rd = c/(2*bw);
lim = ((N_rx-1)*d)^2/lambda;
y_rel = y0/lim;
fprintf('L/rd = %.2f\n', L/rd)
fprintf('Near-far limit = %.2f\n', y_rel)

tic 
image_full = PSF_polar('N_rx', N_rx, 'y0', y0, 'bw',bw,'c', c,'t_p', t_p,'d',d,'fs',fs, 'fc', fc, 'Nu', Nu); % [theta x range]
toc

%%

% Full image 
clims = [-50 0];
f = figure('Position',[360 198 764 420]);
image_full_db = db(abs(image_full.data./max(image_full.data(:))));
imagesc(image_full.u, image_full.range, image_full_db.', clims)
hold on

L_lmb = (N_rx-1)*d/lambda;
lim = ((N_rx-1)*d)^2/lambda;
y_rel = y0/lim;
limit_lmb = ((N_rx-1)*d/lambda)^2;
rd = c/(2*bw);

cb = colorbar(); 
cb.Label.VerticalAlignment = "bottom";
cb.LineWidth = 1;
ylabel(cb,'Power [dB]','FontSize',11,'Rotation',270, 'Interpreter', 'latex')
title('Beamformed image')
% subtitle(sprintf('B/fc = $%.2f$, $%.f$ elements, $d/\\lambda = %.1f$,  $L = %.f\\lambda$, $L^2/\\lambda = %.f\\lambda$, $\\frac{y}{L^2/\\lambda} = %.2f$, F\\# = $%.1f$', bw/fc, N_rx, d/lambda, L_lmb, limit_lmb, y_rel, fnum),'Interpreter','latex')
% subtitle(sprintf('B/fc = $%.2f$, $%.f$ elements, $d/\\lambda = %.1f$, $L = %.f\\lambda$, $L/\\Delta r = %.1f$, $L^2/\\lambda = %.f\\lambda$, $\\frac{y}{L^2/\\lambda} = %.2f$, F\\# = $%.1f$', bw/fc, N_rx, d/lambda, L_lmb, L/rd,  limit_lmb, y_rel, fnum),'Interpreter','latex')
% subtitle(sprintf('B/fc = $%.2f$, $%.f$ elements, $d/\\lambda = %.1f$, $L = %.f\\lambda$, $L/\\Delta r = %.1f$, $\\frac{y}{L^2/\\lambda} = %.2f$, F\\# = $%.1f$', bw/fc, N_rx, d/lambda, L_lmb, L/rd, y_rel, fnum),'Interpreter','latex')
subtitle(sprintf('B/fc = $%.2f$, $%.f$ elements, $d/\\lambda = %.1f$, $L = %.f\\lambda$, $L/\\Delta r = %.1f$, $L^2/\\lambda = %.f\\lambda$, $\\frac{y}{L^2/\\lambda} = %.2f$', bw/fc, N_rx, d/lambda, L_lmb, L/rd,  limit_lmb, y_rel),'Interpreter','latex')

% Plot range spread lines
% hold on
% spread = abs(image_full.u/4*L);
% bw_lim = (abs(image_full.u) < (c/bw)*2/L);
% spread(bw_lim) = c/(2*bw);
% plot(image_full.u, spread+y0, 'r', 'LineWidth', 1)
% plot(image_full.u, -spread+y0, 'r', 'LineWidth', 1)

% k=1;
% xline(k*(c/(fc+bw/2))/d, 'r', 'LineWidth',1, 'LineStyle','-')
% xline(k*(c/(fc-bw/2))/d, 'r', 'LineWidth',1, 'LineStyle','-')
% fnum = y0/L;
% plot(image_full.u, y0*image_full.u/(fnum*2)+y0, 'r', 'LineWidth', 1)
% plot(image_full.u, -y0*image_full.u/(fnum*2)+y0, 'r', 'LineWidth', 1)

xlabel('$u$', 'Interpreter', 'latex')
ylabel('Range [m]')
set(gca,'LineWidth', 1, 'FontName', 'Serif', 'Fontsize', 11)

% figname = sprintf('./Code/figs/grating-lobe-analysis/4/imagefull_rb%.2f_y%.1f_L%.1f_fc%.f_d%.1f.pdf', bw/fc, y0, L, fc, d/lambda);
% exportgraphics(f, figname, 'Resolution', 300)

% figname = sprintf('./Code/figs/grating-lobe-analysis/5/imagefull_rb%.2f_f%.1f_L%.1f_fc%.f_d%.1f.png', bw/fc, fnum, L, fc, d/lambda);
% exportgraphics(f, figname, 'Resolution', 300)
 
% figname = sprintf('./Code/figs/grating-lobe-analysis/u-dependence/imagefull_rb%.2f_f%.1f_L%.1f_fc%.f_d%.1f.png', bw/fc, fnum, L, fc, d/lambda);
% exportgraphics(f, figname, 'Resolution', 300)
% 
% figname = sprintf('./Code/figs/grating-lobe-analysis/u-dependence/imagefull_rb%.2f_f%.1f_L%.1f_fc%.f_d%.1f.pdf', bw/fc, fnum, L, fc, d/lambda);
% exportgraphics(f, figname, 'Resolution', 300)

%%
% PSF 
psf = max(image_full.data, [], 2);

L_lmb = (N_rx-1)*d/lambda;
lim = ((N_rx-1)*d)^2/lambda;
y_rel = y0/lim;
limit_lmb = ((N_rx-1)*d/lambda)^2;
nr_globes = round(d/lambda)-1;

f = figure('Position', [305 241 843 327]);
plot(image_full.u, db(abs(psf./max(psf))), 'LineWidth', 1);
xlabel('$u$', 'Interpreter', 'latex')
ylabel('Power [dB]')
xlim([image_full.u(1) image_full.u(end)])
title('Point spread function')
% subtitle(sprintf('B/fc = $%.2f$, $%.f$ elements, $d/\\lambda = %.1f$,  $L = %.1f\\lambda$, $L^2/\\lambda = %.f\\lambda$, $\\frac{y}{L^2/\\lambda} = %.2f$,  F\\# = $%.1f$', bw/fc, N_rx, d/lambda, L_lmb, limit_lmb, y_rel, y0/L),'Interpreter','latex')
subtitle(sprintf('B/fc = $%.2f$, $%.f$ elements, $d/\\lambda = %.1f$, $L = %.f\\lambda$, $L/\\Delta r = %.1f$, $L^2/\\lambda = %.f\\lambda$, $\\frac{y}{L^2/\\lambda} = %.2f$', bw/fc, N_rx, d/lambda, L_lmb, L/rd,  limit_lmb, y_rel),'Interpreter','latex')

grid on
set(gca,'LineWidth',1, 'FontName', 'Serif', 'Fontsize', 11)
% for i=-nr_globes:1:nr_globes
%     xline(i*lambda/d, 'LineWidth',1, 'LineStyle','--')
% end
% 
% k=1;
% xline(k*(c/(fc+bw/2))/d, 'r', 'LineWidth',1, 'LineStyle','-')
% xline(k*(c/(fc-bw/2))/d, 'r', 'LineWidth',1, 'LineStyle','-')

% figname = sprintf('./Code/figs/grating-lobe-analysis/psf_rb%.2f_y%.f_L%.1f_fc%.f_d%.1f.png', bw/fc, y0, L, fc, d/lambda);
% exportgraphics(f, figname, 'Resolution', 300)

% figname = sprintf('./Code/figs/grating-lobe-analysis/5/psf_rb%.2f_f%.1f_L%.1f_fc%.f_d%.1f.png', bw/fc, fnum, L, fc, d/lambda);
% exportgraphics(f, figname, 'Resolution', 300)

% figname = sprintf('./Code/figs/grating-lobe-analysis/u-dependence/psf_rb%.2f_y%.f_L%.1f_fc%.f_d%.1f.png', bw/fc, y0, L, fc, d/lambda);
% exportgraphics(f, figname, 'Resolution', 300)
% 
% figname = sprintf('./Code/figs/grating-lobe-analysis/u-dependence/psf_rb%.2f_f%.1f_L%.1f_fc%.f_d%.1f.png', bw/fc, fnum, L, fc, d/lambda);
% exportgraphics(f, figname, 'Resolution', 300)

%% 1)
% Far field, undersampled, narrowband, vary freq. 

y0 = 100;
bw = 1e2; fc = 1e4; 
c = 1500; t_p = 10e-2; N_rx = 8;
d = 3*c/fc; Nu = 1000;

lambda1 = c/fc;
fc2 = 2/3*fc;
fcs = [fc fc2];

bw2 = 0.01*fc2;
bws = [bw bw2];

psfs = zeros(length(bws),Nu);
for i = 1:length(fcs)
    image_full = PSF_polar('N_rx', N_rx, 'y0', y0, 'bw',bws(i),'c', c,'t_p', t_p,'d',d,'fs',200e3, 'fc', fcs(i), 'Nu', Nu); % [theta x range]
    psf = max(image_full.data, [], 2);
    psfs(i,:) = psf;
end

%%

f = figure('Position', [360 166 705 450]);
ax1 = subplot(2,1,1);
plot(image_full.u, db(abs(psfs(1,:)./max(psfs(1,:)))), 'LineWidth', 1);
set(gca,'Xticklabel',[])
title('Point spread function')
subtitle(sprintf('B/fc $= %.2f$, $%.f$ elements, $d = %g\\lambda_1$, $L = %.f\\lambda_1$, $L^2/\\lambda_1 = %.f\\lambda_1$, $\\frac{y}{L^2/\\lambda_1} = %.1f$', bw/fc, N_rx,d/lambda1, L_lmb,  limit_lmb, y_rel), 'FontSize', 11, 'FontName', 'Serif', 'Interpreter','latex')
% title(sprintf('$\\lambda_1/d = %.1f$', (c/fcs(1))/d),'Interpreter','latex')
ylabel('Power [dB]')
legend(sprintf('$\\lambda_1/d = %.1f$', (c/fcs(1))/d),'Interpreter','latex')
grid on
for i=[-2 -1 1 2]
    xline(i*(c/fc)/d, 'LineWidth',1, 'LineStyle','--', 'HandleVisibility','off')
end

ax2 = subplot(2,1,2);
plot(image_full.u, db(abs(psfs(2,:)./max(psfs(2,:)))), 'LineWidth', 1);
% title(sprintf('$\\lambda_2/d = %.1f$', (c/fcs(2))/d),'Interpreter','latex')
legend(sprintf('$\\lambda_2/d = %.1f$', (c/fcs(2))/d),'Interpreter','latex')
xlabel('$u$', 'Interpreter', 'latex')
ylabel('Power [dB]')
grid on
for i=[-1 1]
    xline(i*(c/fc2)/d, 'LineWidth',1, 'LineStyle','--', 'HandleVisibility','off')
end

L_lmb = (N_rx-1)*d/lambda1;
y_rel = y0/(image_full.lim);
limit_lmb = ((N_rx-1)*d/lambda1)^2;
linkaxes([ax1,ax2],'xy')
set([ax1,ax2], 'LineWidth', 1, 'FontSize', 11, 'FontName', 'Serif')
% sgtitle(sprintf('\\textbf{Point spread function}\n\\textbf{B/fc $\\mathbf{= %.2f}$, $\\mathbf{%.f}$ elements, $\\mathbf{d = %g\\lambda_1, L = %.f\\lambda_1, L^2/\\lambda_1 = %.f\\lambda_1, \\frac{y}{L^2/\\lambda_1} = %.1f}$}', bw/fc, N_rx,d/lambda1, L_lmb,  limit_lmb, y_rel), 'FontSize', 12, 'FontName', 'Serif', 'Interpreter','latex')

newcolors = [0/255 138/255 222/255;
             255/255 31/255 91/255;
             255/255 198/255 30/255;
             156/255 156/255 156/255;];
         
colororder(newcolors)
 
figname = sprintf('./Code/figs/grating-lobe-analysis/1/far-field-freq_y%.f_L%.1f_rb%.2f.pdf', y0, N_rx*d-d, bw/fc);
exportgraphics(f, figname, 'Resolution', 300)
figname = sprintf('./Code/figs/grating-lobe-analysis/1/far-field-freq_y%.f_L%.1f_rb%.2f.png', y0, N_rx*d-d, bw/fc);
exportgraphics(f, figname, 'Resolution', 300)

%% 2)
% Far field, undersampled, narrowband, fixed freq, vary length

y0 = 200;
bw = 1e2; fc = 1e4; 
c = 1500; t_p = 10e-2;
lambda = c/fc; d = 2*lambda; 
N_rx1 = 8;
N_rx2 = 16;
Ns = [N_rx1 N_rx2];
Nu = 1000;

psfs = zeros(length(N_rx1),1000);
for N=1:length(Ns)
    image_full = PSF_polar('N_rx', Ns(N), 'y0', y0, 'bw',bw,'c', c,'t_p', t_p,'d',d,'fs',200e3, 'fc', fc, 'Nu', Nu); % [theta x range]
    psf = max(image_full.data, [], 2);
    psfs(N,:) = psf;
end

%%

L1_lmb = (N_rx1-1)*d/lambda;
lim = ((N_rx1-1)*d)^2/lambda;
y_rel = y0/lim;
limit_lmb = ((N_rx1-1)*d/lambda)^2;

f = figure('Position', [360 166 705 450]);
ax1 = subplot(2,1,1);
plot(image_full.u, db(abs(psfs(1,:)./max(psfs(1,:)))), 'LineWidth', 1);
set(gca,'Xticklabel',[])
ylabel('Power [dB]')
yline(-13, '-.', '$-13$ dB $\,\,\,$', 'Interpreter', 'latex', 'LineWidth', 1)
title(sprintf('$%.f$ elements, $L = %.f\\lambda$, $L^2/\\lambda = %.f\\lambda$, $\\frac{y}{L^2/\\lambda} = %.1f$', N_rx1, L1_lmb, limit_lmb, y_rel), 'Interpreter','latex')
grid on

for i=[-1 0 1]
    xline(i*(c/fc)/d + (c/fc)/(N_rx1*d), 'LineWidth',1, 'LineStyle','--')
    xline(i*(c/fc)/d - (c/fc)/(N_rx1*d), 'LineWidth',1, 'LineStyle','--')
end

L2_lmb = (N_rx2-1)*d/lambda;
lim = ((N_rx2-1)*d)^2/lambda;
y_rel = y0/lim;
limit_lmb = ((N_rx2-1)*d/lambda)^2;

ax2 = subplot(2,1,2);
plot(image_full.u, db(abs(psfs(2,:)./max(psfs(2,:)))), 'LineWidth', 1);
title(sprintf('$%.f$ elements, $L = %.f\\lambda$, $L^2/\\lambda = %.f\\lambda$, $\\frac{y}{L^2/\\lambda} = %.1f$', N_rx2, L2_lmb, limit_lmb, y_rel), 'Interpreter','latex')
xlabel('$u$', 'Interpreter', 'latex')
ylabel('Power [dB]')
yline(-13, '-.', '$-13$ dB $\,\,\,$', 'Interpreter', 'latex', 'LineWidth', 1)
grid on
for i=[-1 0 1]
    xline(i*(c/fc)/d + (c/fc)/(N_rx2*d), 'LineWidth',1, 'LineStyle','--')
    xline(i*(c/fc)/d - (c/fc)/(N_rx2*d), 'LineWidth',1, 'LineStyle','--')
end

linkaxes([ax1,ax2],'xy')
set([ax1,ax2], 'LineWidth', 1, 'FontSize', 11,'FontName', 'Serif')
sgtitle(sprintf('Point spread function\nB/fc $= %.2f$, $d/\\lambda = %.1f$', bw/fc, d/lambda), 'FontSize', 13, 'FontName', 'Serif', 'Interpreter','latex')

newcolors = [0/255 138/255 222/255;
             255/255 31/255 91/255;
             255/255 198/255 30/255;
             156/255 156/255 156/255;];

colororder(newcolors)

% figname = sprintf('./Code/figs/grating-lobe-analysis/2/far-field-length_y%.f_rb%.2f__N1%.f_N2%.f.png', y0, bw/fc, N_rx1, N_rx2);
% exportgraphics(f, figname, 'Resolution', 300)
% figname = sprintf('./Code/figs/grating-lobe-analysis/2/far-field-length_y%.f_rb%.2f__N1%.f_N2%.f.pdf', y0, bw/fc, N_rx1, N_rx2);
% exportgraphics(f, figname, 'Resolution', 300)

%% 3.0)
% Far field, short array, B/fc = 0.3 example

y0 = 100;
bw = 5e3;
fc = 1e4; 
c = 1500; t_p = 10e-2; N_rx = 8;
lambda = c/fc; d = 3*lambda; 
L = (N_rx-1)*d;
Nu = 500;
nr_globes = 3;

image_full_bws = PSF_polar('N_rx', N_rx, 'y0', y0, 'bw',bw,'c', c,'t_p', t_p,'d',d,'fs',200e3, 'fc', fc, 'Nu', Nu); % [theta x range]
psf = max(image_full_bws.data, [], 2);

%%
% PSF plot

f = figure('Position', [305 241 843 327]);
plot(image_full_bws.u, db(abs(psf./max(psf))), 'LineWidth', 1);
for k=[-nr_globes:-1 1:nr_globes]
    if nr_globes == 3
        % For plotting overlapping grating lobes
        xline(k*lambda/d, 'LineWidth', 1, 'LineStyle','--')
        xline(k*(c/(fc+bw/2))/d, 'LineWidth', 1, 'LineStyle', '-', 'Color', newcolors(abs(k)+1,:))
        xline(k*(c/(fc-bw/2))/d, 'LineWidth', 1, 'LineStyle', '-', 'Color', newcolors(abs(k)+1,:))
    else 
        % For plotting not overlapping grating lobes
        xline(k*lambda/d, 'LineWidth',1, 'LineStyle','--')
        xline(k*(c/(fc+bw/2))/d, 'LineWidth',1, 'LineStyle','-')
        xline(k*(c/(fc-bw/2))/d, 'LineWidth',1, 'LineStyle','-')
    end
end

% Mark main lobe 
xline(lambda/(N_rx*d), 'LineWidth',1, 'LineStyle','-')
xline(-lambda/(N_rx*d), 'LineWidth',1, 'LineStyle','-')

xlabel('$u$', 'Interpreter', 'latex')
ylabel('Power [dB]')
xlim([-1 1])
title('Point spread function')

L_lmb = (N_rx-1)*d/lambda;
lim = ((N_rx-1)*d)^2/lambda;
y_rel = y0/lim;
limit_lmb = ((N_rx-1)*d/lambda)^2;

subtitle(sprintf('B/fc = $%.2f$, $%.f$ elements, $d/\\lambda = %.1f$,  $L = %.f\\lambda$, $L^2/\\lambda = %.f\\lambda$, $\\frac{y}{L^2/\\lambda} = %.1f$', bw/fc, N_rx, d/lambda, L_lmb, limit_lmb, y_rel),'Interpreter','latex')
grid on
set(gca,'LineWidth',1, 'FontSize', 11, 'FontName', 'Serif')

lgd = {};
         
colororder(newcolors)

figname = sprintf('./Code/figs/grating-lobe-analysis/3/psf-BR%.2f.png', bw/fc);
exportgraphics(f, figname, 'Resolution', 300)
figname = sprintf('./Code/figs/grating-lobe-analysis/3/psf-BR%.2f.pdf', bw/fc);
exportgraphics(f, figname, 'Resolution', 300)

%% 3.1)
% Far field, short array, relative BW vary
y0 = 50;
% bws = linspace(1e2,6e3,60);  % for 3.5k ++ blir estimatet feil.
bws = [2e3 6e3 8e3 1e4];
fc = 1e4; 
c = 1500; t_p = 10e-2; N_rx = 8;
L = N_rx*d-d;
Nu = 1000;
lambda = c/fc; d = 2*lambda; 
nr_globes = round(d/lambda) - 1; % Nr. of grating lobes in visible region

psfs = zeros(length(bws),Nu);
grels = zeros(length(bws),nr_globes);

for i = 1:length(bws)
    fprintf('B = %.f Hz \n', bws(i));
    image_full_bws = PSF_polar('N_rx', N_rx, 'y0', y0, 'bw',bws(i),'c', c,'t_p', t_p,'d',d,'fs',200e3, 'fc', fc, 'Nu', Nu); % [theta x range]
    psfs(i,:) = max(image_full_bws.data, [], 2);
    grels(i, :) = REL_LEVEL(psfs(i,:), lambda, d, c, fc,bws(i), image_full_bws.u, N_rx*d-d, N_rx); % Relative grating lobe level
end

%% 3.1) Plot of PSFs
% Vary relative bandwidths

f = figure('Position', [360 191 817 318]);
hold on

for i = 1:length(bws)
    plot(image_full_bws.u, db(abs(psfs(i,:)./max(psfs(i,:)))), 'LineWidth', 1)
    lgd{i} = sprintf('B/fc = %.2f', bws(i)/fc);
end

L_lmb = (N_rx-1)*d/lambda;
lim = ((N_rx-1)*d)^2/lambda;
y_rel = y0/lim;
limit_lmb = ((N_rx-1)*d/lambda)^2;

% Mark main lobes
xline(lambda/(N_rx*d), 'LineWidth', 1, 'LineStyle', '--')
xline(-lambda/(N_rx*d), 'LineWidth', 1, 'LineStyle', '--')

title(sprintf('Point spread functions'))
subtitle(sprintf('$%.f$ elements, $d/\\lambda = %.1f$,  $L = %.f\\lambda$, $L^2/\\lambda = %.f\\lambda$, $\\frac{y}{L^2/\\lambda} = %.1f$',  N_rx, d/lambda, L_lmb, limit_lmb, y_rel),'Interpreter','latex')
set(gca,'LineWidth', 1, 'Fontsize', 11, 'FontName', 'Serif')
grid on
xlim([-0.2530 0.2530])
xlabel('$u$', 'Interpreter', 'latex')
ylabel('Power [dB]')
legend(lgd, 'Location', 'bestoutside')

colororder(newcolors)
lgd = {};

% figname = sprintf('./Code/figs/grating-lobe-analysis/3/psfs_rb_yrel%.1f_L%.1f_d%.1f.png', y_rel, L, d/lambda);
% exportgraphics(f, figname, 'Resolution', 300)
% figname = sprintf('./Code/figs/grating-lobe-analysis/3/psfs_rb_yrel%.1f_L%.1f_d%.1f.pdf', y_rel, L, d/lambda);
% exportgraphics(f, figname, 'Resolution', 300)

%% 3.1)
% Plot relative grating lobe level as a function of relative bandwidth

f = figure('Position',[360 212 653 298]);
hold on
plot(bws/fc, grels(:,1), 'LineWidth', 1)
legend('p $= 1$', 'Interpreter', 'latex')

title('Relative grating lobe level')
L_lmb = (N_rx-1)*d/lambda;
lim = ((N_rx-1)*d)^2/lambda;
y_rel = y0/lim;
limit_lmb = ((N_rx-1)*d/lambda)^2;

subtitle(sprintf('$%.f$ elements, $d/\\lambda = %.1f$,  $L = %.f\\lambda$, $L^2/\\lambda = %.f\\lambda$, $\\frac{y}{L^2/\\lambda} = %.1f$', N_rx, d/lambda, L_lmb, limit_lmb, y_rel),'Interpreter','latex')
xlabel('B/fc')
ylabel('Power [dB]')
set(gca,'LineWidth', 1, 'Fontsize', 11, 'FontName', 'Serif')
grid on

% figname = sprintf('./Code/figs/grating-lobe-analysis/3/grel_rb_d%.1f_y%.f_L%.1f.png', d/lambda, y0, N_rx*d-d);
% exportgraphics(f, figname, 'Resolution', 300)
% figname = sprintf('./Code/figs/grating-lobe-analysis/3/grel_rb_d%.1f_y%.f_L%.1f.pdf', d/lambda, y0, N_rx*d-d);
% exportgraphics(f, figname, 'Resolution', 300)

%% 3.1) + Angelsen
% Plot relative grating lobe level as a function of relative bandwidth

f = figure('Position',[360 212 653 298]);
hold on
plot(bws/fc, grels(:,1), 'LineWidth', 1)
angels= db(1.33/(N_rx*1)*fc./bws);
plot(bws/fc, angels, 'LineWidth', 1)

legend('p $= 1$','Angelsen, p $= 1$', 'Interpreter', 'latex')
title('Relative grating lobe level')
L_lmb = (N_rx-1)*d/lambda;
lim = ((N_rx-1)*d)^2/lambda;
y_rel = y0/lim;
limit_lmb = ((N_rx-1)*d/lambda)^2;

subtitle(sprintf('$%.f$ elements, $d/\\lambda = %.1f$,  $L = %.f\\lambda$, $L^2/\\lambda = %.f\\lambda$, $\\frac{y}{L^2/\\lambda} = %.1f$', N_rx, d/lambda, L_lmb, limit_lmb, y_rel),'Interpreter','latex')
xlabel('B/fc')
ylabel('Power [dB]')
set(gca,'LineWidth', 1, 'Fontsize', 11, 'FontName', 'Serif')
grid on

newcolors = [0/255 138/255 222/255;
             255/255 31/255 91/255;
             255/255 198/255 30/255;
             156/255 156/255 156/255;];
         
colororder(newcolors)

figname = sprintf('./Code/figs/grating-lobe-analysis/3/grel_rb_angelsen_d%.1f_y%.f_L%.1f.png', d/lambda, y0, N_rx*d-d);
exportgraphics(f, figname, 'Resolution', 300)
figname = sprintf('./Code/figs/grating-lobe-analysis/3/grel_rb_angelsen_d%.1f_y%.f_L%.1f.pdf', d/lambda, y0, N_rx*d-d);
exportgraphics(f, figname, 'Resolution', 300)



%% 3.2)
% Far field, different array lengths, relative BW vary

y0 = 60;
bws = linspace(1e2,5e3,50);  % for 3.5k ++ blir estimatet feil.
fc = 1e4; 
c = 1500; t_p = 10e-2; N_rx = 4;
L = N_rx*d-d;
Nu = 500;
lambda = c/fc; d = 2.5*lambda; 
nr_globes = round(d/lambda) - 1; % Nr. of grating lobes in visible region

N_rx2 = 8;
Ns = [N_rx N_rx2];

psfs = zeros(length(Ns), length(bws),Nu);
grels = zeros(length(Ns), length(bws),nr_globes);

for n=1:length(Ns)
    fprintf('N = %.f\n', Ns(n));
    for i = 1:length(bws)
        fprintf('B = %.f Hz \n', bws(i));
        image_full = PSF_polar('N_rx', Ns(n), 'y0', y0, 'bw',bws(i),'c', c,'t_p', t_p,'d',d,'fs',200e3, 'fc', fc, 'Nu', Nu); % [theta x range]
        psfs(n, i,:) = max(image_full.data, [], 2);
        grels(n, i, :) = REL_LEVEL(psfs(n,i,:), lambda, d, c, fc,bws(i), image_full.u, Ns(n)*d-d, Ns(n)); % Relative grating lobe level
    end
end

%% 3.2)
% plot

f = figure('Position',[360 212 653 298]);
hold on
for n=1:length(Ns)
    for g=1:nr_globes
        plot(bws(1:40)/fc, grels(n,1:40,g), 'LineWidth', 1)
    end 
end

% plot(bws(1:40)/fc, db(1.33/(Ns(1)*1)*fc./bws(1:40))) % Angelsen
lgd1 = sprintf('$N_1 = %.1f$, p=1', N_rx);
lgd12 = sprintf('$N_1 = %.1f$, p=2', N_rx);
lgd2 = sprintf('$N_2 = %.1f$, p=1', N_rx2);
lgd22 = sprintf('$N_2 = %.1f$, p=2', N_rx2);

legend(lgd1, lgd12, lgd2, lgd22, 'Interpreter', 'latex', 'Location', 'bestoutside')
L_lmb1 = (Ns(1)-1)*d/lambda;
L_lmb2 = (Ns(2)-1)*d/lambda;
lim = ((Ns(2)-1)*d)^2/lambda;
y_rel = y0/lim;
limit_lmb = ((Ns(2)-1)*d/lambda)^2;

title('Relative grating lobe level')
subtitle(sprintf('$d/\\lambda = %.1f$,  $L_1 = %.f\\lambda$,  $L_2 = %.f\\lambda$, $\\frac{y}{L_{2}^2/\\lambda} = %.1f$', d/lambda, L_lmb1, L_lmb2, y_rel),'Interpreter','latex')
xlabel('B/fc')
ylabel('Power [dB]')
set(gca,'LineWidth', 1, 'Fontsize', 11, 'FontName', 'Serif')
grid on

% figname = sprintf('./Code/figs/grating-lobe-analysis/3/grel_rb_d%.1f_y%.f_N1-%.f_N2-%.f.pdf', d/lambda, y0, Ns(1), Ns(2));
% exportgraphics(f, figname, 'Resolution', 300)

%% 3.3)
% Far field, different element spacings, relative BW vary, long run time

y0 = 100;
bws = linspace(1e2,2.8e3,28);  % for 3.5k ++ blir estimatet feil.
% bws = linspace(1e2,6e3,10);  % for 3.5k ++ blir estimatet feil.
fc = 1e4; 
c = 1500; t_p = 10e-2; N_rx = 8;
Nu = 500;
lambda = c/fc; d = 2.5*lambda; d2 = 3.5*lambda;
% d = 1.5*lambda; d2 = 5*lambda;
L2 = N_rx*d2-d2;
nr_globes1 = round(d/lambda) - 1; % Nr. of grating lobes in visible region
nr_globes2 = round(d2/lambda) - 1; % Nr. of grating lobes in visible region

ds = [d d2];
%%
psfs = zeros(length(ds), length(bws),Nu);
grels = zeros(length(bws),nr_globes1);
grels2 = zeros(length(bws),nr_globes2);

for j=1:length(ds)
    fprintf('d = %.3f\n', ds(j));
    for i = 1:length(bws)
        fprintf('B = %.f Hz \n', bws(i));
        image_full = PSF_polar('N_rx', N_rx, 'y0', y0, 'bw',bws(i),'c', c,'t_p', t_p,'d',ds(j),'fs',200e3, 'fc', fc, 'Nu', Nu); % [theta x range]
        psfs(j, i,:) = max(image_full.data, [], 2);
        if j == 1
            grels(i, :) = REL_LEVEL(psfs(j,i,:), lambda, ds(j), c, fc,bws(i), image_full.u, N_rx*ds(j)-ds(j), N_rx); % Relative grating lobe level
        else
            grels2(i, :) = REL_LEVEL(psfs(j,i,:), lambda, ds(j), c, fc,bws(i), image_full.u, N_rx*ds(j)-ds(j), N_rx); % Relative grating lobe level
        end
    end
end

%% 3.3) plot
% Different relative grating lobe levels
f = figure('Position',[360 191 730 319]);
hold on

j = 1;
for i=1:nr_globes1
    plot(bws/fc, grels(:,i), 'LineWidth', 1)
%     lgd{j} = sprintf('$d_1/\\lambda = %.1f$, p $= %.f$', d/lambda, i);
    lgd{j} = sprintf('$u_1 = \\lambda/d_1 = %.1f$, p $= %.f$', lambda/d, i);
    j = j +1;
end

for i=1:nr_globes2
    plot(bws/fc, grels2(:,i), 'LineWidth', 1)
%     lgd{j} = sprintf('$d_2/\\lambda = %.1f$, p $= %.f$', d2/lambda, i);
%     lgd{j} = sprintf('$d_2/\\lambda = %.1f$, p $= %.f$', d2/lambda, i);
    lgd{j} = sprintf('$u_2 = \\lambda/d_2 = %.1f$, p $= %.f$', lambda/d2, i);
    j = j + 1;
end

legend(lgd, 'Interpreter', 'latex', 'Location','best')

title('Relative grating lobe level')
L_lmbmax = (N_rx-1)*d2/lambda;
lim = ((N_rx-1)*d2)^2/lambda;
y_rel = y0/lim;
limit_lmb_max = ((N_rx-1)*d2/lambda)^2;

xlim([bws(1)/fc bws(end)/fc])

subtitle(sprintf('$%.f$ elements, $L_{max} = %.f\\lambda$, $L_{max}^2/\\lambda = %.f\\lambda$, $\\frac{y}{L_{max}^2/\\lambda} = %.1f$', N_rx, L_lmbmax, limit_lmb_max, y_rel),'Interpreter','latex')
xlabel('B/fc')
ylabel('Power [dB]')
set(gca,'LineWidth', 1, 'Fontsize', 11, 'FontName', 'Serif')
grid on

newcolors = [255/255 31/255 91/255;
             0/255 138/255 222/255;
             0/255 205/255 109/255;
             255/255 198/255 30/255;
             242/255 133/252 34/255];
         
colororder(newcolors)

lgd = {};
% figname = sprintf('./Code/figs/grating-lobe-analysis/3/grel_rb_dmin%.1f_dmax%.1f_y%.f.png', d/lambda, d2/lambda, y0);
% exportgraphics(f, figname, 'Resolution', 300)
% figname = sprintf('./Code/figs/grating-lobe-analysis/3/grel_rb_dmin%.1f_dmax%.1f_y%.f.pdf', d/lambda, d2/lambda, y0);
% exportgraphics(f, figname, 'Resolution', 300)

%% 3.3) plot with angelsen
f = figure('Position',[360 182 836 327]);
hold on

j = 1;
for i=1:nr_globes1
    plot(bws/fc, grels(:,i), 'LineWidth', 1)
    lgd{j} = sprintf('$d_1/\\lambda = %.1f$, p $= %.f$', d/lambda, i);
    j = j +1;
end

for i=1:nr_globes2
    plot(bws/fc, grels2(:,i), 'LineWidth', 1)
    lgd{j} = sprintf('$d_2/\\lambda = %.1f$, p $= %.f$', d2/lambda, i);
    j = j + 1;
end

% Angelsen
for i=1:nr_globes2
    plot(bws/fc, db(1.33/(N_rx*i)*fc./bws), 'LineWidth', 1) 
    lgd{j} = sprintf('Angelsen, p $= %.f$', i);
    j = j + 1;
end

legend(lgd, 'Interpreter', 'latex', 'Location','bestoutside')
title('Relative grating lobe level')
L_lmbmax = (N_rx-1)*d2/lambda;
lim = ((N_rx-1)*d2)^2/lambda;
y_rel = y0/lim;
limit_lmb_max = ((N_rx-1)*d2/lambda)^2;

subtitle(sprintf('$%.f$ elements, $L_{max} = %.f\\lambda$, $L_{max}^2/\\lambda = %.f\\lambda$, $\\frac{y}{L_{max}^2/\\lambda} = %.1f$', N_rx, L_lmbmax, limit_lmb_max, y_rel),'Interpreter','latex')
xlabel('B/fc')
ylabel('Power [dB]')
set(gca,'LineWidth', 1, 'Fontsize', 11, 'FontName', 'Serif')
grid on

xlim([bws(1)/fc bws(end)/fc])

lgd = {};
% figname = sprintf('./Code/figs/grating-lobe-analysis/3/grel_rb_angelsen_dmin%.1f_dmax%.1f_y%.f.png', d/lambda, d2/lambda, y0);
% exportgraphics(f, figname, 'Resolution', 300)
% figname = sprintf('./Code/figs/grating-lobe-analysis/3/grel_rb_angelsen_dmin%.1f_dmax%.1f_y%.f.pdf', d/lambda, d2/lambda, y0);
% exportgraphics(f, figname, 'Resolution', 300)

%% 3.3) plot with multiplying
f = figure('Position',[360 191 730 319]);
hold on

j = 1;
for i=1:nr_globes1
    plot(bws*i/fc, grels(:,i), 'LineWidth', 1)
    if i == 1
       lgd{j} = sprintf('$d_1/\\lambda = %.1f$, p $= %.f$', d/lambda, i);
    else 
       lgd{j} = sprintf('$d_1/\\lambda = %.1f$, p $= %.f$\nmult. by %.1f', d/lambda, i, i);
    end
    j = j +1;
end

for i=1:nr_globes2
    plot(bws*i/fc, grels2(:,i), 'LineWidth', 1)
    if i == 1
       lgd{j} = sprintf('$d_2/\\lambda = %.1f$, p $= %.f$', d2/lambda, i);
    else
       lgd{j} = sprintf('$d_2/\\lambda = %.1f$,p $= %.f$\nmult. by %.1f', d2/lambda, i, i);
    end 
    j = j + 1;
end

legend(lgd, 'Interpreter', 'latex', 'Location','best')
title('Relative grating lobe level')

L_lmbmax = (N_rx-1)*d2/lambda;
lim = ((N_rx-1)*d2)^2/lambda;
y_rel = y0/lim;
limit_lmb_max = ((N_rx-1)*d2/lambda)^2;

subtitle(sprintf('$%.f$ elements, $L_{max} = %.f\\lambda$, $L_{max}^2/\\lambda = %.f\\lambda$, $\\frac{y}{L_{max}^2/\\lambda} = %.1f$', N_rx, L_lmbmax, limit_lmb_max, y_rel),'Interpreter','latex')
xlabel('B/fc')
ylabel('Power [dB]')
set(gca,'LineWidth', 1, 'Fontsize', 11, 'FontName', 'Serif')
grid on

xlim([bws(1)/fc bws(end)*3/fc])


newcolors = [255/255 31/255 91/255;
             0/255 138/255 222/255;
             0/255 205/255 109/255;
             255/255 198/255 30/255;
             242/255 133/252 34/255];
         
colororder(newcolors)


lgd = {};
figname = sprintf('./Code/figs/grating-lobe-analysis/3/grel_rb_mult_dmin%.1f_dmax%.1f_y%.f.png', d/lambda, d2/lambda, y0);
exportgraphics(f, figname, 'Resolution', 300)
figname = sprintf('./Code/figs/grating-lobe-analysis/3/grel_rb_mult_dmin%.1f_dmax%.1f_y%.f.pdf', d/lambda, d2/lambda, y0);
exportgraphics(f, figname, 'Resolution', 300)
%% 3.4) 
% Near field, vary relative bandwidth

% bws = linspace(1e2,10e3,5);
bws = [2e3 6e3 8e3 1e4];
fc = 1e4; 
c = 1500; t_p = 10e-2; N_rx = 8;
L = N_rx*d-d;
Nu = 1000;
fnum = 10;
y0 = fnum*L;
lambda = c/fc; d = 2*lambda; 
nr_globes = round(d/lambda) - 1; % Nr. of grating lobes in visible region

psfs = zeros(length(bws),Nu);
grels = zeros(length(bws),nr_globes);

for i = 1:length(bws)
    fprintf('B = %.f Hz \n', bws(i));
    image_full_bws = PSF_polar('N_rx', N_rx, 'y0', y0, 'bw',bws(i),'c', c,'t_p', t_p,'d',d,'fs',200e3, 'fc', fc, 'Nu', Nu); % [theta x range]
    psfs(i,:) = max(image_full_bws.data, [], 2);
    grels(i, :) = REL_LEVEL(psfs(i,:), lambda, d, c, fc,bws(i), image_full_bws.u, N_rx*d-d, N_rx); % Relative grating lobe level
end
%% 3.4) Plot
% Near field, vary relative bandwidth

f = figure('Position', [360 191 817 318]);
hold on

for i = 1:length(bws)
    plot(image_full_bws.u, db(abs(psfs(i,:)./max(psfs(i,:)))), 'LineWidth', 1)
    lgd{i} = sprintf('B/fc = %.2f', bws(i)/fc);
end

L_lmb = (N_rx-1)*d/lambda;
lim = ((N_rx-1)*d)^2/lambda;
y_rel = y0/lim;
limit_lmb = ((N_rx-1)*d/lambda)^2;

% Mark main lobes
xline(lambda/(N_rx*d), 'LineWidth', 1, 'LineStyle', '--')
xline(-lambda/(N_rx*d), 'LineWidth', 1, 'LineStyle', '--')

title(sprintf('Point spread functions\n$%.f$ elements, $d/\\lambda = %.1f$, $L = %.f\\lambda$, $L^2/\\lambda = %.f\\lambda$, $\\frac{y}{L^2/\\lambda} = %.1f$, F\\# $= %.1f$', N_rx, d/lambda, L_lmb, limit_lmb, y_rel, fnum),'Interpreter','latex')
set(gca,'LineWidth', 1, 'Fontsize', 11, 'FontName', 'Serif')
grid on
xlim([-0.2530 0.2451])
xlabel('$u$', 'Interpreter', 'latex')
ylabel('Power [dB]')
legend(lgd, 'Location', 'bestoutside')

figname = sprintf('./Code/figs/grating-lobe-analysis/3/psfs_rb_yrel%.1f_L%.1f_d%.1f.png', y_rel, L, d/lambda);
exportgraphics(f, figname, 'Resolution', 300)
figname = sprintf('./Code/figs/grating-lobe-analysis/3/psfs_rb_yrel%.1f_L%.1f_d%.1f.pdf', y_rel, L, d/lambda);
exportgraphics(f, figname, 'Resolution', 300)

lgd = {};
%% 4.1) 
% Far field, relative bandwidth, array length vary

y0 = 250;
fc = 1e4; bw = 3e3;
c = 1500; t_p = 10e-2; 
lambda = c/fc; d = 2.5*lambda; 
Nu = 1000;
nr_globes = round(d/lambda) - 1; % Nr. of grating lobes in visible region
drange = c/(2*bw);
Nr = length((y0-10):drange/10:(y0+10));
Ns = [4 6 8 10 12 14 16];
% Ns = 16;
Ns = [5 10 15];

psfs = zeros(length(Ns), Nu);
grels = zeros(length(Ns),nr_globes);
imgs = zeros(Nu, Nr, length(Ns));

for j=1:length(Ns)
    fprintf('N = %.f elements\n', Ns(j))
    image_full = PSF_polar("bw",bw,'fc', fc, 'N_rx', Ns(j), 'Nu', Nu, 't_p', t_p, 'fs',200e3, 'y0', y0,'d', d);  % [theta x range]
    imgs(:,:,j) = image_full.data;
    psfs(j,:) = max(image_full.data, [], 2);
    grels(j, :) = REL_LEVEL(psfs(j,:), lambda, d, c, fc, bw, image_full.u, Ns(j)*d-d, Ns(j)); % Relative grating lobe level
end


%% 4.1) Plot the PSFs
% as a function of elements

f = figure('Position', [305 241 843 327]);

hold on
for j=1:length(Ns)
    plot(image_full.u, db(abs(psfs(j,:)./max(psfs(j,:)))), 'LineWidth', 1);
    lgd{j} = sprintf('N: %.f', Ns(j));
end

L_max = Ns(end)*d-d;
Lmax_lmb = L_max/lambda;
max_limit_lmb = ((Ns(end)-1)*d/lambda)^2;
y_rel = y0/(L_max^2/lambda);
rd = c/(2*bw);

legend(lgd, 'Location', 'bestoutside')
xlabel('$u$', 'Interpreter','latex')
ylabel('Power [dB]')
title('Point spread function')
subtitle(sprintf('B/fc = $%.2f$, $\\Delta r = %.2f$ m, $d/\\lambda = %.1f$, $L_{max} = %.f\\lambda$, $L_{max}^2/\\lambda = %.f\\lambda$, $\\frac{y}{L_{max}^2/\\lambda} = %.1f$', bw/fc, rd, d/lambda, Lmax_lmb, max_limit_lmb, y_rel),'Interpreter','latex')
grid on
set(gca,'LineWidth',1, 'FontName', 'Serif', 'Fontsize', 11)

% figname = sprintf('./Code/figs/grating-lobe-analysis/psfs_length_%.1f-%.1f.png', Ns(1)*d-d, Ns(end)*d-d);
% exportgraphics(f, figname, 'Resolution', 300)
% figname = sprintf('./Code/figs/grating-lobe-analysis/psfs_length_%.1f-%.1f.pdf', Ns(1)*d-d, Ns(end)*d-d);
% exportgraphics(f, figname, 'Resolution', 300)

lgd = {};

%% 4.1) Range spread: PSFs vs PW beampattern
% as a function of elements

grels_pw_bp = zeros(length(Ns),nr_globes);
pw_beampatterns = zeros(length(Ns), 8192);
for j=1:length(Ns)
    fprintf('N = %.f elements\n', Ns(j))
    [total_BP, angles, u] = PW_BP("bw",bw, 'c',c, 'N',Ns(j), 'fc', fc, 'd', d, 'NFFT', size(pw_beampatterns, 2), 'umax', 0.99);
    pw_beampatterns(j, :) = total_BP;
    grels_pw_bp(j,:) = REL_LEVEL(total_BP, lambda, d, c, fc, bw, u, Ns(j)*d-d, Ns(j));
end

%% 4.1) Range spread: PSFs vs PW beampattern
% Plot PSFs
f = figure('Position', [607 1007 892 536]);

L_max = Ns(end)*d-d;
Lmax_lmb = L_max/lambda;
max_limit_lmb = ((Ns(end)-1)*d/lambda)^2;
y_rel = y0/(L_max^2/lambda);
rd = c/(2*bw);
sgtitle(sprintf('Range spread effect\nB/fc = $%.2f$, $\\Delta r = %.2f$ m, $L_{max}/\\Delta r = %.2f$, $d/\\lambda = %.1f$', bw/fc, rd, (Ns(end)*d-d)/rd, d/lambda),'Interpreter','latex', 'Fontsize', 12)

subplot(2,1,1)
hold on
k = 1;
for j=1:length(Ns)
    plot(image_full.u, db(abs(psfs(j,:)./max(psfs(j,:)))), 'LineWidth', 1);
    lgd{k} = sprintf('N: %.f', Ns(j));
    k = k+1;
end
legend(lgd, 'Location', 'bestoutside')
xlabel('$u$', 'Interpreter','latex')
ylabel('Power [dB]')
title('Point spread function')
subtitle(sprintf('$\\frac{y}{L_{max}^2/\\lambda} = %.1f$', y_rel),'Interpreter','latex')

grid on
set(gca,'LineWidth',1, 'FontName', 'Serif', 'Fontsize', 11)
lgd = {};

subplot(2,1,2)
hold on
k = 1;
for j=1:length(Ns)
    plot(u, db(abs(pw_beampatterns(j,:)./max(pw_beampatterns(j,:)))), 'LineWidth', 1);
    lgd{k} = sprintf('N: %.f', Ns(j));
    k= k+1;
end
legend(lgd, 'Location', 'bestoutside')
xlabel('$u$', 'Interpreter','latex')
ylabel('Power [dB]')
title('Pulsed wave beampattern')
grid on
ylim([-30 0])
set(gca,'LineWidth',1, 'FontName', 'Serif', 'Fontsize', 11)

lgd = {};

% figname = sprintf('./Code/figs/grating-lobe-analysis/4/pwbeampattern-bwpsf.png');
% exportgraphics(f, figname, 'Resolution', 300)
% figname = sprintf('./Code/figs/grating-lobe-analysis/4/pwbeampattern-bwpsf.pdf');
% exportgraphics(f, figname, 'Resolution', 300)

%% 4.1) Range spread: PSFs vs PW beampattern
% Relative grating lobe level, plot

angelsen1 = db(1.33*fc./(Ns*bw));
angelsen2 = db(1.33*fc./(2*Ns*bw));

f = figure('Position', [751 907 843 511]);
sgtitle(sprintf('Relative grating lobe level\nB/fc = $%.2f$, $\\Delta r = %.2f$ m, $L_{max}/\\Delta r = %.2f$, $d/\\lambda = %.1f$', bw/fc, rd, (Ns(end)*d-d)/rd, d/lambda),'Interpreter','latex','Fontsize', 12)


% Relative grating lobe level
subplot(2,1,1)
plot(Ns, grels_pw_bp(:,1), 'LineWidth', 1)
hold on
plot(Ns, grels(:,1), 'LineWidth', 1)
plot(Ns, angelsen1, 'LineWidth', 1)
legend('PW beampattern', 'PSF', 'Angelsen', 'Interpreter', 'latex')
title('$|p| = 1$', 'Interpreter', 'latex')
xlabel('Number of elements')
ylabel('Power [dB]')
grid on
set(gca,'LineWidth',1, 'FontName', 'Serif', 'Fontsize', 11)


subplot(2,1,2)
plot(Ns, grels_pw_bp(:,2), 'LineWidth', 1)
hold on
plot(Ns, grels(:,2), 'LineWidth', 1)
plot(Ns, angelsen2, 'LineWidth', 1)
legend('PW beampattern', 'PSF', 'Angelsen', 'Interpreter', 'latex')
title('$|p| = 2$', 'Interpreter', 'latex')
xlabel('Number of elements')
ylabel('Power [dB]')
grid on
set(gca,'LineWidth',1, 'FontName', 'Serif', 'Fontsize', 11)

% figname = sprintf('./Code/figs/grating-lobe-analysis/4/pwbeampattern-bwpsf-grel.png');
% exportgraphics(f, figname, 'Resolution', 300)
% figname = sprintf('./Code/figs/grating-lobe-analysis/4/pwbeampattern-bwpsf-grel.pdf');
% exportgraphics(f, figname, 'Resolution', 300)
%% 4.1)
% Plot relative grating lobe level, (at max)
f = figure('Position', [305 241 843 327]);

% Relative grating lobe level
plot(Ns, grels(:,1), 'LineWidth', 1)
hold on
plot(Ns, grels(:,2), 'LineWidth', 1)
legend('p $= 1$', 'p $= 2$', 'Interpreter', 'latex')
title('Relative grating lobe level')
xlabel('Number of elements')
ylabel('Power [dB]')
grid on
subtitle(sprintf('B/fc = $%.2f$, $\\Delta r = %.2f$ m, $d/\\lambda = %.1f$, $L_{max} = %.f\\lambda$, $L_{max}^2/\\lambda = %.f\\lambda$, $\\frac{y}{L_{max}^2/\\lambda} = %.1f$', bw/fc, rd, d/lambda, Lmax_lmb, max_limit_lmb, y_rel),'Interpreter','latex')

set(gca,'LineWidth',1, 'FontName', 'Serif', 'Fontsize', 11)
% figname = sprintf('./Code/figs/grating-lobe-analysis/4/grel_length_%.1f-%.1f.png', Ns(1)*d-d, Ns(end)*d-d);
% exportgraphics(f, figname, 'Resolution', 300)

% Relative grating lobe level mult. by grating lobe number
f = figure('Position', [305 241 843 327]);

plot(Ns, grels(:,1), 'LineWidth', 1)
hold on
plot(Ns*2, grels(:,2), 'LineWidth', 1)
legend('p $= 1$', sprintf('p $= 2$\nmult. by 2'), 'Interpreter', 'latex')
title('Relative grating lobe level')
xlabel('Number of elements')
ylabel('Power [dB]')
grid on
subtitle(sprintf('B/fc = $%.2f$, $\\Delta r = %.2f$ m, $d/\\lambda = %.1f$, $L_{max} = %.f\\lambda$, $L_{max}^2/\\lambda = %.f\\lambda$, $\\frac{y}{L_{max}^2/\\lambda} = %.1f$', bw/fc, rd, d/lambda, Lmax_lmb, max_limit_lmb, y_rel),'Interpreter','latex')

set(gca,'LineWidth',1, 'FontName', 'Serif', 'Fontsize', 11)
% figname = sprintf('./Code/figs/grating-lobe-analysis/4/grel_length_mult_%.1f-%.1f.png', Ns(1)*d-d, Ns(end)*d-d);
% exportgraphics(f, figname, 'Resolution', 300)

% Relative grating lobe level compare with Angelsen
f = figure('Position', [305 241 843 327]);

angelsen1 = db(1.33*fc./(Ns*bw));
angelsen2 = db(1.33*fc./(2*Ns*bw));
plot(Ns, grels(:,1), 'LineWidth', 1)
hold on
plot(Ns, grels(:,2), 'LineWidth', 1)
plot(Ns, angelsen1, 'LineWidth', 1)
plot(Ns, angelsen2, 'LineWidth', 1)
legend('p $= 1$', 'p $= 2$', 'Angelsen, p=1', 'Angelsen, p=2', 'Interpreter', 'latex')
title('Relative grating lobe level')
xlabel('Number of elements')
ylabel('Power [dB]')
grid on
subtitle(sprintf('B/fc = $%.2f$, $\\Delta r = %.2f$ m, $d/\\lambda = %.1f$, $L_{max} = %.f\\lambda$, $L_{max}^2/\\lambda = %.f\\lambda$, $\\frac{y}{L_{max}^2/\\lambda} = %.1f$', bw/fc, rd, d/lambda, Lmax_lmb, max_limit_lmb, y_rel),'Interpreter','latex')

set(gca,'LineWidth',1, 'FontName', 'Serif', 'Fontsize', 11)
% figname = sprintf('./Code/figs/grating-lobe-analysis/4/grel_length_angelsen_%.1f-%.1f.png', Ns(1)*d-d, Ns(end)*d-d);
% exportgraphics(f, figname, 'Resolution', 300)

%% 4.1) 
% Save full beamformed image for different array lengths

for j=1:length(Ns)
    f = figure('Position',[360 198 764 420]);
    clims = [-50 0];
    image_full_db = db(abs(imgs(:,:,j)./max(imgs(:,:,j), [], 'all')));
    imagesc(image_full.u, image_full.range, image_full_db.', clims)
    
    L = Ns(j)*d-d;
    L_lmb = L/lambda;
    limit_lmb = ((Ns(j)-1)*d/lambda)^2;
    y_rel = y0/(L^2/lambda);
    rd = c/(2*bw);

    cb = colorbar(); 
    cb.Label.VerticalAlignment = "bottom";
    cb.LineWidth = 1;
    ylabel(cb,'Power [dB]','FontSize',11,'Rotation',270, 'FontName', 'Serif')
    title('Beamformed image')
    subtitle(sprintf('B/fc = $%.2f$, %.f elements, $d/\\lambda = %.1f$, $L = %.f\\lambda$, $L/\\Delta r = %.2f$, $L^2/\\lambda = %.f\\lambda$, $\\frac{y}{L^2/\\lambda} = %.1f$', bw/fc, Ns(j),d/lambda, L_lmb, L/rd,  limit_lmb, y_rel),'Interpreter','latex')

    xlabel('$u$', 'Interpreter', 'latex')
    ylabel('Range [m]')
    set(gca,'LineWidth', 1, 'FontName', 'Serif', 'Fontsize', 11)

%     figname = sprintf('./Code/figs/grating-lobe-analysis/4/imagefull_rb%.2f_y%.f_L%.1f_fc%.f_d%.1f.png', bw/fc, y0, Ns(j)*d-d, fc, d/lambda);
%     exportgraphics(f, figname, 'Resolution', 300)
% 
%     figname = sprintf('./Code/figs/grating-lobe-analysis/4/imagefull_rb%.2f_y%.f_L%.1f_fc%.f_d%.1f.pdf', bw/fc, y0, Ns(j)*d-d, fc, d/lambda);
%     exportgraphics(f, figname, 'Resolution', 300)
end

%% 4.1)
% Calculate rel 2d globe level

save = 0; % Save image ?
draw = 0; % Plot + draw trapezoid ?
avg = 1; % Calculate relative average power ?

rel_2dglobe = zeros(length(Ns), nr_globes+1);
for j=1:length(Ns)
    rel_2dglobe(j,:) = rel_blob_power(imgs(:,:,j),Ns(j),y0,fc,c,bw,lambda,image_full.range,image_full.u,d,save, draw,avg);
end

%% 4.1)
% Plot relative grating lobe average power inside boxes
f = figure('Position', [305 241 843 327]);

L = Ns(end)*d-d;
Lmax_lmb = L/lambda;
max_limit_lmb = ((Ns(j)-1)*d/lambda)^2;
y_rel = y0/(L^2/lambda);
rd = c/(2*bw);

plot(Ns, rel_2dglobe(:,1), 'LineWidth', 1)
hold on
plot(Ns, rel_2dglobe(:,2), 'LineWidth', 1)
legend('p $= 1$', 'p $= 2$', 'Interpreter', 'latex')
title('Relative average grating lobe level')
xlabel('Number of elements')
ylabel('Power [dB]')
grid on
subtitle(sprintf('B/fc = $%.2f$, $\\Delta r = %.2f$ m, $d/\\lambda = %.1f$, $L_{max} = %.f\\lambda$, $L_{max}^2/\\lambda = %.f\\lambda$, $\\frac{y}{L_{max}^2/\\lambda} = %.1f$', bw/fc, rd, d/lambda, Lmax_lmb, max_limit_lmb, y_rel),'Interpreter','latex')

set(gca,'LineWidth',1, 'FontName', 'Serif', 'Fontsize', 11)

% figname = sprintf('./Code/figs/grating-lobe-analysis/4/grel_avg_%.1f-%.1f.png', Ns(1)*d-d, Ns(end)*d-d);
% exportgraphics(f, figname, 'Resolution', 300)
% figname = sprintf('./Code/figs/grating-lobe-analysis/4/grel_avg_%.1f-%.1f.pdf', Ns(1)*d-d, Ns(end)*d-d);
% exportgraphics(f, figname, 'Resolution', 300)

f = figure('Position', [305 241 843 327]);

plot(Ns, rel_2dglobe(:,1), 'LineWidth', 1)
hold on
plot(Ns*2, rel_2dglobe(:,2), 'LineWidth', 1)
legend('p $= 1$', sprintf('p $= 2$\nmult. by 2'), 'Interpreter', 'latex')
title('Relative average grating lobe level')
xlabel('Number of elements')
ylabel('Power [dB]')
grid on
subtitle(sprintf('B/fc = $%.2f$, $\\Delta r = %.2f$ m, $d/\\lambda = %.1f$, $L_{max} = %.f\\lambda$, $L_{max}^2/\\lambda = %.f\\lambda$, $\\frac{y}{L_{max}^2/\\lambda} = %.1f$', bw/fc, rd, d/lambda, Lmax_lmb, max_limit_lmb, y_rel),'Interpreter','latex')

set(gca,'LineWidth',1, 'FontName', 'Serif', 'Fontsize', 11)

% figname = sprintf('./Code/figs/grating-lobe-analysis/4/grel_avg_mult_%.1f-%.1f.png', Ns(1)*d-d, Ns(end)*d-d);
% exportgraphics(f, figname, 'Resolution', 300)
% figname = sprintf('./Code/figs/grating-lobe-analysis/4/grel_avg_mult_%.1f-%.1f.pdf', Ns(1)*d-d, Ns(end)*d-d);
% exportgraphics(f, figname, 'Resolution', 300)

f = figure('Position', [305 241 843 327]);

angelsen1 = db(1.33*fc./(Ns*bw));
angelsen2 = db(1.33*fc./(2*Ns*bw));
plot(Ns, rel_2dglobe(:,1), 'LineWidth', 1)
hold on
plot(Ns, rel_2dglobe(:,2), 'LineWidth', 1)
plot(Ns, angelsen1, 'LineWidth', 1)
plot(Ns, angelsen2, 'LineWidth', 1)
legend('p $= 1$', 'p $= 2$', 'Angelsen, p = 1', 'Angelsen, p = 2', 'Interpreter', 'latex')
title('Relative average grating lobe level')
xlabel('Number of elements')
ylabel('Power [dB]')
grid on
subtitle(sprintf('B/fc = $%.2f$, $\\Delta r = %.2f$ m, $d/\\lambda = %.1f$, $L_{max} = %.f\\lambda$, $L_{max}^2/\\lambda = %.f\\lambda$, $\\frac{y}{L_{max}^2/\\lambda} = %.1f$', bw/fc, rd, d/lambda, Lmax_lmb, max_limit_lmb, y_rel),'Interpreter','latex')

set(gca,'LineWidth',1, 'FontName', 'Serif', 'Fontsize', 11)
% 
% figname = sprintf('./Code/figs/grating-lobe-analysis/4/grel_avg_angelsen_%.1f-%.1f.png', Ns(1)*d-d, Ns(end)*d-d);
% exportgraphics(f, figname, 'Resolution', 300)
% figname = sprintf('./Code/figs/grating-lobe-analysis/4/grel_avg_angelsen_%.1f-%.1f.pdf', Ns(1)*d-d, Ns(end)*d-d);
% exportgraphics(f, figname, 'Resolution', 300)

%% 4.1)
% Plot sum of power inside a blob for mainlobe+grating lobes

save = 0; % Save image ?
draw = 0; % Plot + draw trapezoid ?
avg = 0; % Calculate relative average power ?

rel_2dglobe = zeros(length(Ns), nr_globes+1);
for j=1:length(Ns)
    rel_2dglobe(j,:) = rel_blob_power(imgs(:,:,j),Ns(j),y0,fc,c,bw,lambda,image_full.range,image_full.u,d,save, draw,avg);
end
%%
f = figure('Position', [305 267 744 301]);

L = Ns(end)*d-d;
Lmax_lmb = L/lambda;
max_limit_lmb = ((Ns(j)-1)*d/lambda)^2;
y_rel = y0/(L^2/lambda);

plot(Ns, rel_2dglobe(:,1), 'LineWidth', 1)
hold on
plot(Ns, rel_2dglobe(:,2), 'LineWidth', 1)
plot(Ns, rel_2dglobe(:,3), 'LineWidth', 1)
legend('Mainlobe', 'p $= 1$', 'p $= 2$', 'Interpreter', 'latex', 'Location', 'best')
title('Sum of power inside trapezoids')
xlabel('Number of elements')
ylabel('Power')
grid on
subtitle(sprintf('B/fc = $%.2f$, $d/\\lambda = %.1f$, $L_{max} = %.f\\lambda$, $L_{max}^2/\\lambda = %.f\\lambda$, $\\frac{y}{L_{max}^2/\\lambda} = %.1f$', bw/fc, d/lambda, Lmax_lmb, max_limit_lmb, y_rel),'Interpreter','latex')

set(gca,'LineWidth',1, 'FontName', 'Serif', 'Fontsize', 11)

colororder(newcolors)

figname = sprintf('./Code/figs/grating-lobe-analysis/4/gsum_trapezoid_%.1f-%.1f.png', Ns(1)*d-d, Ns(end)*d-d);
exportgraphics(f, figname, 'Resolution', 300)
figname = sprintf('./Code/figs/grating-lobe-analysis/4/gsum_trapezoid_%.1f-%.1f.pdf', Ns(1)*d-d, Ns(end)*d-d);
exportgraphics(f, figname, 'Resolution', 300)


%% 5) Near field 

fnums = linspace(1,10,100);
bw = 1e2; fc = 1e4; 
c = 1500; t_p = 10e-2; N_rx = 8;
Nu = 500;
lambda = c/fc; d = 4*lambda; 
L = d*(N_rx-1);
y = fnums.*L;
nr_globes = round(d/lambda) - 1; % Nr. of grating lobes 

%%
psfs = zeros(length(y), Nu);
grels = zeros(length(y),nr_globes);

for j=1:length(y)
    fprintf('fnum = %.1f \n', fnums(j))
    image_full = PSF_polar("bw",bw,'fc', fc, 'N_rx', N_rx, 'Nu', Nu, 't_p', t_p, 'fs',200e3, 'y0', y(j),'d', d);  % [theta x range]
    psfs(j,:) = max(image_full.data, [], 2);
    grels(j, :) = REL_LEVEL(psfs(j,:), lambda, d, c, fc, bw, image_full.u, d*(N_rx-1), N_rx); % Relative grating lobe level
end

%% 5.1) 
% Relative grating lobe level

% f = figure('Position', [305 148.3333 710 392.6667]);
f = figure('Position', [304 269 715 299]);
hold on
for i=1:nr_globes
    plot(fnums, grels(:,i), 'LineWidth', 1)
end

legend('p $= 1$', 'p $= 2$', 'p $= 3$', 'Location', 'bestoutside')
xlabel('F-number')
ylabel('Power [dB]')
grid on

L_lmb = L/lambda;
limit_lmb = ((N_rx-1)*d/lambda)^2;
title('Relative grating lobe level')
titletxt = sprintf('B/fc = $%.2f$, $%.f$ elements, $d/\\lambda = %.1f$,  $L = %.f\\lambda$, $L^2/\\lambda = %.f\\lambda$', bw/fc, N_rx, d/lambda, L_lmb, limit_lmb);
subtitle(titletxt,'Interpreter','latex',  'FontSize', 12);

set(gca,'LineWidth', 1, 'Fontsize', 11, 'FontName', 'Serif')
grid on

% figname = sprintf('./Code/figs/grating-lobe-analysis/5/fnum_y_d%.1f_fnum%.1f.png', d/lambda, fnum(end));
% exportgraphics(f, figname, 'Resolution', 300)
% figname = sprintf('./Code/figs/grating-lobe-analysis/5/fnum_y_d%.1f_fnum%.1f.pdf', d/lambda, fnum(end));
% exportgraphics(f, figname, 'Resolution', 300)

lgd = {};

%% 5.1)
% Try to plot the different PSFs for different ranges 

f = figure('Position', [360 198 919 420])
k=1;
for i=[1 2 3 10]
    plot(image_full.u, db(abs(psfs(i, :)/max(psfs(i, :), [], 'all'))), 'LineWidth', 1)
    hold on
    lgd{k} = sprintf('F-number = %.1f', fnums(i));
    k = k+1;
end

set(gca,'LineWidth', 1, 'Fontsize', 11, 'FontName', 'Serif')
grid on
legend(lgd,'Location', 'bestoutside')
title('Point spread function')
subtitle(sprintf('$B/f_c = %.2f$, $L^2/\\lambda = %.1f$ m, $%.f$ elements, L = $%.1f$ m, $d/\\lambda = %.1f$', bw/fc, image_full.lim, N_rx, L, d/lambda),'Interpreter','latex');

lgd = {};

% figname = sprintf('./Code/figs/grating-lobe-analysis/5/psfs_y_d%.1f_rb%.2f_L%.1f.png', d/lambda, bw/fc, N_rx*d-d);
% exportgraphics(f, figname, 'Resolution', 300)
%% 5.2) Near field 
% Two different number of elements

fnums = linspace(1,5,50);
bw = 1e2; fc = 1e4; 
c = 1500; t_p = 10e-2; 
N_rx1 = 4; N_rx2 = 8;
Nu = 500;
lambda = c/fc; d = 4*lambda; 
L1 = d*(N_rx1-1);
L2 = d*(N_rx2-1);
y1 = fnums.*L1;
y2 = fnums.*L2;
nr_globes = round(d/lambda) - 1; % Nr. of grating lobes 
Ns = [N_rx1 N_rx2];
ys = [y1; y2];

psfs = zeros(length(Ns), length(y), Nu);
grels = zeros(length(Ns), length(y),nr_globes);
% imgs = zeros(Nu, Nr, length(y));

%%
for n=1:length(Ns)
    for j=1:length(fnums)
        fprintf('fnum = %.1f\n', fnums(j))
        image_full = PSF_polar("bw",bw,'fc', fc, 'N_rx', Ns(n), 'Nu', Nu, 't_p', t_p, 'fs',200e3, 'y0', ys(n,j),'d', d);  % [theta x range]
        psfs(n,j,:) = max(image_full.data, [], 2);
        grels(n,j, :) = REL_LEVEL(psfs(n,j,:), lambda, d, c, fc, bw, image_full.u, d*(Ns(n)-1), Ns(n)); % Relative grating lobe level
    end
end

%% 5.2)
% Relative grating lobe level

f = figure('Position', [304 218 751 349]);
hold on

j = 1;
for n=1:length(Ns)
    for i=1:nr_globes
        plot(fnums, grels(n,:,i), 'LineWidth', 1)
        lgd{j} = sprintf('N$_%.f = %.f$, p $= %.f$', n, Ns(n), i);
        j=j+1;
    end
end

legend(lgd,'Interpreter', 'latex', 'Location', 'bestoutside')
xlabel('F-number')
ylabel('Power [dB]')
grid on

title('Relative grating lobe level')
titletxt = sprintf('B/fc = $%.2f$, $d/\\lambda = %.1f$, $L_1 = %.1f\\lambda$, $L_2 = %.1f\\lambda$', bw/fc, d/lambda, d*(N_rx1-1)/lambda, d*(N_rx2-1)/lambda);
subtitle(titletxt,'Interpreter','latex',  'FontSize', 12);
set(gca,'LineWidth', 1, 'Fontsize', 11, 'FontName', 'Serif')
grid on
 
% figname = sprintf('./Code/figs/grating-lobe-analysis/5/fnum_y_d%.1f_N1-%.f_N2-%.f_fnum%.1f.png', d/lambda, N_rx1, N_rx2, fnum(end));
% exportgraphics(f, figname, 'Resolution', 300)
% figname = sprintf('./Code/figs/grating-lobe-analysis/5/fnum_y_d%.1f_N1-%.f_N2-%.f_fnum%.1f.pdf', d/lambda, N_rx1, N_rx2, fnum(end));
% exportgraphics(f, figname, 'Resolution', 300)

lgd = {};

%% 6) U-dependence
% Without range spread
% Constant F# number

bw = 1e2; fc = 1e4; 
c = 1500; t_p = 10e-2; 
Nu = 1000;
lambda = c/fc; 
u = linspace(0.1, 0.9, 30);
d = lambda./u;
Lc = 28*lambda;
N_rx = round(Lc./d+1);

L = N_rx.*d-d;
% fnums = 10;
fnums = [2.5 5 10];
y0 = fnums.'*L; % [fnums x lengths]
rd = c/(2*bw);

%%
psfs = zeros(length(d), Nu);
grel1 = zeros(length(fnums),length(d)); % Relative globe level of first grating lobe

for f=1:length(fnums)
    for j=1:length(d)
        fprintf('d/lambda = %.1f, u = %.2f\n', d(j)/lambda, lambda/d(j))
        image_full = PSF_polar('N_rx', N_rx(j), 'y0', y0(f,j), 'bw', bw,'c', c,'t_p', t_p,'d',d(j),'fs',200e3, 'fc', fc, 'Nu', Nu);
        psfs(j,:) = max(image_full.data, [], 2);
        grel = REL_LEVEL(psfs(j,:), lambda, d(j), c, fc, bw, image_full.u, L(j), N_rx(j)); % Relative grating lobe level
        grel1(f,j) = grel(1);
    end
end

%%

f = figure('Position', [0 1 734 319]);
hold on
for i=1:length(fnums)
    plot(u, grel1(i, :), '-', 'Linewidth', 1)
    lgd{i} = sprintf('F\\# $= %.1f$',fnums(i));
end

xlabel('$u$', 'Interpreter', 'latex')
ylabel('Power [dB]')
legend(lgd, 'Location', 'best','Interpreter', 'latex')
Lmax = mean(L); %(N_rx-1)*d(1);
rd = c/(2*bw);

title('Relative grating lobe level')
titletxt = sprintf('B/fc = $%.2f$, $L = %.1f\\lambda$, $L/\\Delta r = %.2f$', bw/fc, Lmax/lambda, Lmax/rd);
subtitle(titletxt,'Interpreter','latex',  'FontSize', 12);
set(gca,'LineWidth', 1, 'Fontsize', 11, 'FontName', 'Serif')
grid on
lgd = {};

newcolors = [0/255 138/255 222/255;
            255/255 31/255 91/255;
             0/255 205/255 109/255;
             255/255 198/255 30/255;
             242/255 133/252 34/255];
colororder(newcolors)

figname = sprintf('./Code/figs/grating-lobe-analysis/u-dependence/grel_fnum1%.1f_fnum2%.1f_rb%.2f_L%.1f.png', fnums(1),fnums(end), bw/fc, Lmax);
exportgraphics(f, figname, 'Resolution', 300)
figname = sprintf('./Code/figs/grating-lobe-analysis/u-dependence/grel_fnum1%.1f_fnum2%.1f_rb%.2f_L%.1f.pdf', fnums(1),fnums(end), bw/fc, Lmax);
exportgraphics(f, figname, 'Resolution', 300)
%% Vary the array length and see if any different
% Near field
% Without range spread 

bw = 1e2; fc = 1e4; 
c = 1500; t_p = 10e-2;
Nu = 1000;
lambda = c/fc; 
u = linspace(0.1, 0.9, 30);
d = lambda./u;
L1w = 28*lambda; % Set length 
L2w = 14*lambda; % Set length
N_rx1 = round(28*lambda./d+1);
N_rx2 = round(14*lambda./d+1);

L1 = N_rx1.*d-d;
L2 = N_rx2.*d-d;
fnums = 2.5;
% fnums = [2.5 5 10 15];
y1 = fnums.'*L1; % [fnums x lengths]
y2 = fnums.'*L2;
rd = c/(2*bw);

%%
psfs = zeros(length(d), Nu);
grel1 = zeros(length(fnums),length(d)); 
grel2 = zeros(length(fnums),length(d)); 

for f=1:length(fnums)
    for j=1:length(d)
        fprintf('d/lambda = %.1f, u = %.2f\n', d(j)/lambda, lambda/d(j))
        image_full1 = PSF_polar('N_rx', N_rx1(j), 'y0', y1(f,j), 'bw', bw,'c', c,'t_p', t_p,'d',d(j),'fs',200e3, 'fc', fc, 'Nu', Nu);
        psfs(j,:) = max(image_full1.data, [], 2);
        grel = REL_LEVEL(psfs(j,:), lambda, d(j), c, fc, bw, image_full.u, L1(j), N_rx1(j)); % Relative grating lobe level
        grel1(f,j) = grel(1);

        image_full2 = PSF_polar('N_rx', N_rx2(j), 'y0', y2(f,j), 'bw', bw,'c', c,'t_p', t_p,'d',d(j),'fs',200e3, 'fc', fc, 'Nu', Nu);
        psfs(j,:) = max(image_full2.data, [], 2);
        grel = REL_LEVEL(psfs(j,:), lambda, d(j), c, fc, bw, image_full.u, L2(j), N_rx2(j)); % Relative grating lobe level
        grel2(f,j) = grel(1);
    end
end

%%

f = figure('Position', [383 245 754 304]);
% figure;
hold on
for i=1:length(fnums)
    plot(u, grel1(i, :), '-', 'Linewidth', 1)
    plot(u, grel2(i, :), '-', 'Linewidth', 1)
%     lgd{i} = sprintf('L = %.1f,F\\# $= %.1f$',fnums(i));
end

legend(sprintf('$L_1 = %.f\\lambda$', L1w/lambda), sprintf('$L_2 = %.f\\lambda$', L2w/lambda))
xlabel('$u$', 'Interpreter', 'latex')
ylabel('Power [dB]')
legend(lgd, 'Location', 'bestoutside','Interpreter', 'latex')
rd = c/(2*bw);

title('Relative grating lobe level', 'FontSize', 12, 'FontName', 'Serif', 'Interpreter','latex')
titletxt = sprintf('B/fc = $%.2f$, F\\# $= %.1f$, $L_1/\\Delta r = %.1f$', bw/fc, fnums, L1w/rd);
subtitle(titletxt,'Interpreter','latex',  'FontSize', 12);
set(gca,'LineWidth', 1, 'Fontsize', 11, 'FontName', 'Serif')
grid on
lgd = {};

% figname = sprintf('./Code/figs/grating-lobe-analysis/u-dependence/grel_fnum%.1f_L1%.1f_L2%.1f_rb%.2f.png', fnums, L1w, L2w, bw/fc);
% exportgraphics(f, figname, 'Resolution', 300)
% figname = sprintf('./Code/figs/grating-lobe-analysis/u-dependence/grel_fnum%.1f_L1%.1f_L2%.1f_rb%.2f.pdf', fnums, L1w, L2w, bw/fc);
% exportgraphics(f, figname, 'Resolution', 300)

%% 6.1) U-dependence
% With range spread
% Far field

bw = 3e3; fc = 1e4; 
c = 1500; t_p = 10e-2; N_rx = 8;
Nu = 500;
lambda = c/fc; 

% 0.1+1/(1.3043)*(c/(fc-bw/2)-c/(fc+bw/2)); % Fant minste mulige u for høyre
% side av gitterlobe, gitt en båndbredde
% u = linspace(0.1353, 0.9, 30); % U-position to the right of first grating lobe
u = linspace(0.1, 1/(1+1/lambda*(c/(fc-bw/2)-c/(fc+bw/2))/2),30);
% u = linspace(0.35, 1/(1+1/lambda*(c/(fc-bw/2)-c/(fc+bw/2))/2),30);
d = lambda./u; 
Lc = 20*lambda; % Wanted length
% Lc2 = 23*lambda; 
N_rx = round(Lc./d+1);
% N_rx2 = round(Lc2./d+1);

L = N_rx.*d-d;
rd = c/(2*bw);
fprintf('L/rd = %.2f\n', L/rd)
y0 = 80;

drange = c/(2*bw);
Nr = length((y0-10):drange/10:(y0+10));

%%
psfs = zeros(length(d), Nu);
grel1 = zeros(length(d), 1); % Relative globe level of first grating lobe
grel2 = zeros(length(d), 1); 
avg_power = zeros(length(d), 2);
imgs = zeros(length(d), Nu, Nr);

save = 0;
draw = 0;
avg = 1;

for j=1:length(d)
    fprintf('d/lambda = %.1f, u = %.2f\n', d(j)/lambda, lambda/d(j))
    image_full = PSF_polar('N_rx', N_rx(j), 'y0', y0, 'bw', bw,'c', c,'t_p', t_p,'d',d(j),'fs',200e3, 'fc', fc, 'Nu', Nu);
    imgs(j, :, :) = image_full.data;
%     image_full = PSF_polar('N_rx', N_rx, 'y0', y0, 'bw', bw,'c', c,'t_p', t_p,'d',d(j),'fs',200e3, 'fc', fc, 'Nu', Nu);
    psfs(j,:) = max(image_full.data, [], 2);
    grel = REL_LEVEL(psfs(j,:), lambda, d(j), c, fc, bw, image_full.u, (N_rx(j)-1)*d(j), N_rx(j)); % Relative grating lobe level, changing number of elements
    rel_power = rel_blob_power(image_full.data, N_rx(j), y0, fc,c,bw,lambda,image_full.range, image_full.u,d(j), save, draw, avg);
    avg_power(j,:) = rel_power(1:2);
%     grel = REL_LEVEL(psfs(j,:), lambda, d(j), c, fc, bw, image_full.u, L(j), N_rx); % Relative grating lobe level, keeping number of
    % elements constant.
    grel1(j) = grel(1); % Want only the first grating lobe

%     image_full = PSF_polar('N_rx', N_rx2(j), 'y0', y0, 'bw', bw,'c', c,'t_p', t_p,'d',d(j),'fs',200e3, 'fc', fc, 'Nu', Nu);
%     psfs(j,:) = max(image_full.data, [], 2);
%     grel = REL_LEVEL(psfs(j,:), lambda, d(j), c, fc, bw, image_full.u, (N_rx2(j)-1)*d(j), N_rx2(j)); % Relative grating lobe level, changing number of elements
%     grel2(j) = grel(1); % Want only the first grating lobe
end


%%
f = figure('Position', [0 1 734 319]);
hold on

% p = polyfit(u,grel1,2);
% yfit = polyval(p,u);

% lambda_min = c/(fc+bw/2);
p = 1;
% rel_amp = zeros(1,length(N_rx));
% for i=1:length(N_rx)
%     rel_amp(i) = parallel_side_amp(N_rx(i),c,fc,bw,p);
% end
% test = db(rel_amp);

% Calculate rel. average power
rel_avg_power = avg_power(:,2)./avg_power(:,1);

% Calculate areas of trapezoids
% globe_areas = zeros(1, length(N_rx));
[globe_areas, ~, ~, ~] = trapezoid_area(N_rx, bw, fc, c, u, d);
mainlobe_areas = c/bw*lambda./(d.*N_rx)*2;

rel_areas = mainlobe_areas./globe_areas;

% test = db((c/bw)./((N_rx-1)*lambdamin*p/2));

% angels = db(1.33./(N_rx.*1)*fc./bw);
angels = db(1.33*fc*lambda./(d.*u.*N_rx*bw));

plot(u, grel1, '-', 'Linewidth', 1)
% plot(u, grel2, '-', 'Linewidth', 1)
% plot(u, angels, '-', 'Linewidth', 1)
% plot(u, test, '-', 'Linewidth', 1)
plot(u, 10*log10(rel_avg_power), '-', 'Linewidth', 1)
plot(u, 10*log10(rel_areas), '-', 'Linewidth', 1)
% plot(u, yfit, '-', 'Linewidth', 1)
legend('Max', 'Avg. power', 'Area')
% legend('Simulation result', 'Regression fit')
% legend('Simulation result', 'Angelsen', 'Mine')
% legend('Simulation result', 'Mine')

% legend(sprintf('$L_1 = %.f\\lambda$', Lc/lambda), sprintf('$L_2 = %.f\\lambda$', Lc2/lambda), 'Interpreter', 'latex')
xlabel('$u$', 'Interpreter', 'latex')
ylabel('Power [dB]')
rd = c/(2*bw);

title('Relative grating lobe level')
titletxt = sprintf('B/fc = $%.2f$, $L/\\Delta r = %.2f$, $\\frac{y}{L^2/\\lambda} = %.2f$', bw/fc, Lc/rd, y0/(Lc^2/lambda));
% titletxt = sprintf('B/fc = $%.2f$, %.f elements, $L_{min}/\\Delta r = %.2f$, $\\frac{y}{L_{min}^2/\\lambda} = %.2f$', bw/fc, N_rx, min(L)/rd, y0/(Lc^2/lambda));
subtitle(titletxt,'Interpreter','latex',  'FontSize', 12);
set(gca,'LineWidth', 1, 'Fontsize', 11, 'FontName', 'Serif')
grid on
lgd = {};

newcolors = [0/255 138/255 222/255;
            255/255 31/255 91/255;
             0/255 205/255 109/255;
             255/255 198/255 30/255;
             242/255 133/252 34/255];
colororder(newcolors)

% figname = sprintf('./Code/figs/grating-lobe-analysis/u-dependence/grel_far_rb%.2f_L1%.1f_L2%.1f.png', bw/fc, Lc/lambda, Lc2/lambda);
% exportgraphics(f, figname, 'Resolution', 300)
% figname = sprintf('./Code/figs/grating-lobe-analysis/u-dependence/grel_far_rb%.2f_L1%.1f_L2%.1f.pdf', bw/fc, Lc/lambda, Lc2/lambda);
% exportgraphics(f, figname, 'Resolution', 300)

% figname = sprintf('./Code/figs/grating-lobe-analysis/u-dependence/grel_far_rb%.2f_N%.1f.png', bw/fc, N_rx);
% exportgraphics(f, figname, 'Resolution', 300)
% figname = sprintf('./Code/figs/grating-lobe-analysis/u-dependence/grel_far_rb%.2f_N%.1f.pdf', bw/fc, N_rx);
% exportgraphics(f, figname, 'Resolution', 300)

% figname = sprintf('./Code/figs/grating-lobe-analysis/u-dependence/grel_area-vs-pow_rb%.2f.png', bw/fc);
% exportgraphics(f, figname, 'Resolution', 300)
% figname = sprintf('./Code/figs/grating-lobe-analysis/u-dependence/grel_area-vs-pow_rb%.2f.pdf', bw/fc);
% exportgraphics(f, figname, 'Resolution', 300)

%% With and without range spread comparison subplots
% 2d PSFs

clims = [-50 0];
bw = 7e2; fc = 1e4; 
c = 1500; t_p = 10e-2;
Nu = 500;
lambda = c/fc; d = 1.3*lambda; 

N_rx1 = 5;
N_rx2 = 32;
fprintf('image 1\n')
image_full1 = PSF_polar('N_rx', N_rx1, 'y0', y0, 'bw',bw,'c', c,'t_p', t_p,'d',d,'fs',200e3, 'fc', fc, 'Nu', Nu); % [theta x range]
fprintf('image 2\n')
image_full2 = PSF_polar('N_rx', N_rx2, 'y0', y0, 'bw',bw,'c', c,'t_p', t_p,'d',d,'fs',200e3, 'fc', fc, 'Nu', Nu); % [theta x range]

%%
f = figure('Position', [669.6667 926.3333 644.6666 798.7000]);

% Params
L1 = (N_rx1-1)*d;
L_lmb = L1/lambda;
lim = L1^2/lambda;
limit_lmb = lim/lambda;
y_rel = y0/lim;
rd = c/(2*bw);

ax1 = subplot(2,1,1);
image_full_db = db(abs(image_full1.data./max(image_full1.data(:))));
imagesc(image_full1.u, image_full1.range, image_full_db.', clims)

title(sprintf('$%.f$ elements, $L = %.f\\lambda$, $L/\\Delta r = %.1f$, $L^2/\\lambda = %.f\\lambda$, $\\frac{y}{L^2/\\lambda} = %.2f$', N_rx1, L_lmb, L1/rd,  limit_lmb, y_rel),'Interpreter','latex')
xlabel('$u$', 'Interpreter', 'latex')
ylabel('Range [m]')

cb = colorbar(); 
cb.Label.VerticalAlignment = "bottom";
cb.LineWidth = 1;
ylabel(cb,'Power [dB]','FontSize',11,'Rotation',270, 'Interpreter', 'latex')

ax2 = subplot(2,1,2);
image_full_db = db(abs(image_full2.data./max(image_full2.data(:))));
imagesc(image_full2.u, image_full2.range, image_full_db.', clims)

% Plot cross
hold on
plot(image_full.u, image_full.u/4*L2+y0, 'r', 'LineWidth', 1)
plot(image_full.u, -image_full.u/4*L2+y0, 'r', 'LineWidth', 1)

cb = colorbar(); 
cb.Label.VerticalAlignment = "bottom";
cb.LineWidth = 1;
ylabel(cb,'Power [dB]','FontSize',11,'Rotation',270, 'Interpreter', 'latex')

% Params
L2 = (N_rx2-1)*d;
L_lmb = L2/lambda;
lim = L2^2/lambda;
limit_lmb = lim/lambda;
y_rel = y0/lim;

title(sprintf('$%.f$ elements, $L = %.f\\lambda$, $L/\\Delta r = %.1f$, $L^2/\\lambda = %.f\\lambda$, $\\frac{y}{L^2/\\lambda} = %.2f$', N_rx2, L_lmb, L2/rd,  limit_lmb, y_rel),'Interpreter','latex')
xlabel('$u$', 'Interpreter', 'latex')
ylabel('Range [m]')

sgtitle(sprintf('Beamformed images\nB/fc = $%.2f$, $d/\\lambda = %.1f$', bw/fc, d/lambda), 'FontSize', 12, 'FontName', 'Serif', 'Interpreter','latex')
set([ax1,ax2], 'LineWidth', 1, 'FontSize', 11, 'FontName', 'Serif')

% figname = sprintf('./Code/figs/grating-lobe-analysis/4/rangespreadx_rb%.2f_y%.1f_N1-%g_N2-%g_d%.1f.png', bw/fc, y0, N_rx1, N_rx2, d/lambda);
% exportgraphics(f, figname, 'Resolution', 300)
% figname = sprintf('./Code/figs/grating-lobe-analysis/4/rangespreadx_rb%.2f_y%.1f_N1-%g_N2-%g_d%.1f.pdf', bw/fc, y0, N_rx1, N_rx2, d/lambda);
% exportgraphics(f, figname, 'Resolution', 300)

%% U-dependence example subplots 
% 2d PSFs
% Range spread

bw = 3e3; fc = 1e4; 
c = 1500; t_p = 10e-2; N_rx1 = 16; %N_rx1 = 8; %N_rx1 = 16;
Nu = 500;
lambda = c/fc;

u = 0.4;%0.25;%0.4; % U-position of first grating lobe
d1 = lambda/u;
u = 0.8;%0.5;%0.8; % U-position of first grating lobe
d2 = lambda/u;
L = N_rx1*d1-d1;

N_rx2 = round(L/d2+1);

% fnum = 2.5;
% y0 = fnum*L;
y0 = 250;

rd = c/(2*bw);
lim = (L^2/lambda);
y_rel = y0/lim;
fprintf('L/rd = %.2f\n', L/rd)
fprintf('Near-far limit = %.2f\n', y_rel)

L_lmb = L/lambda;
limit_lmb = lim/lambda;
clims = [-50 0];

image_full1 = PSF_polar('N_rx', N_rx1, 'y0', y0, 'bw',bw,'c', c,'t_p', t_p,'d', d1,'fs',200e3, 'fc', fc, 'Nu', Nu); % [theta x range]
image_full2 = PSF_polar('N_rx', N_rx2, 'y0', y0, 'bw',bw,'c', c,'t_p', t_p,'d', d2,'fs',200e3, 'fc', fc, 'Nu', Nu); % [theta x range]

%%
f = figure('Position', [0 0 644.6666 798.7000]);

ax1 = subplot(2,1,1);

image_full_db = db(abs(image_full1.data./max(image_full1.data(:))));
imagesc(image_full1.u, image_full1.range, image_full_db.', clims)

title(sprintf('$%.f$ elements, $d/\\lambda = %.1f$', N_rx1, d1/lambda),'Interpreter','latex')
xlabel('$u$', 'Interpreter', 'latex')
ylabel('Range [m]')

cb = colorbar(); 
cb.Label.VerticalAlignment = "bottom";
cb.LineWidth = 1;
ylabel(cb,'Power [dB]','FontSize',11,'Rotation',270, 'Interpreter', 'latex')

ax2 = subplot(2,1,2);

image_full_db = db(abs(image_full2.data./max(image_full2.data(:))));
imagesc(image_full2.u, image_full2.range, image_full_db.', clims)

title(sprintf('$%.f$ elements, $d/\\lambda = %.1f$', N_rx2, d2/lambda),'Interpreter','latex')
xlabel('$u$', 'Interpreter', 'latex')
ylabel('Range [m]')
sgtitle(sprintf('Beamformed images\n B/fc $= %.2f$, $L = %.f\\lambda$, $L/\\Delta r = %.1f$, $L^2/\\lambda = %.f\\lambda$, $\\frac{y}{L^2/\\lambda} = %.2f$', bw/fc, L_lmb, L/rd,  limit_lmb, y_rel), 'FontSize', 12, 'FontName', 'Serif', 'Interpreter','latex')
% sgtitle(sprintf('Beamformed images\n B/fc $= %.2f$, $L = %.f\\lambda$, $L/\\Delta r = %.1f$, $L^2/\\lambda = %.f\\lambda$, $\\frac{y}{L^2/\\lambda} = %.2f$, $F\\# = %.1f$', bw/fc, L_lmb, L/rd,  limit_lmb, y_rel, fnum), 'FontSize', 12, 'FontName', 'Serif', 'Interpreter','latex')

set([ax1,ax2], 'LineWidth', 1, 'FontSize', 10, 'FontName', 'Serif')

cb = colorbar(); 
cb.Label.VerticalAlignment = "bottom";
cb.LineWidth = 1;
ylabel(cb,'Power [dB]','FontSize',11,'Rotation',270, 'Interpreter', 'latex')

% figname = sprintf('./Code/figs/grating-lobe-analysis/u-dependence/ex_rb%.2f_N1-%g_d1-%.1f_yrel%.1f.png', bw/fc, N_rx1, d1/lambda, y_rel);
% exportgraphics(f, figname, 'Resolution', 300)
% figname = sprintf('./Code/figs/grating-lobe-analysis/u-dependence/ex_rb%.2f_N1-%g_d1-%.1f_yrel%.1f.pdf', bw/fc, N_rx1, d1/lambda, y_rel);
% exportgraphics(f, figname, 'Resolution', 300)

%% U-dependence plots
% PSFs for constant L/rd

bw1 = 1e3; fc = 1e4; 
bw2 = 2e3;
c = 1500; t_p = 10e-2; N_rx = 8; 
Nu = 500;
lambda = c/fc;

lambda_max1 = c./(fc-bw1/2);
lambda_max2 = c./(fc-bw2/2);
d1 = 6*lambda;
d2 = 3*lambda;
L = N_rx*d1-d1;
L2 = N_rx*d2-d2;

y0 = 300;

rd1 = c/(2*bw1);
rd2 = c/(2*bw2);
lim = (L^2/lambda);
y_rel = y0/lim;
fprintf('L/rd = %.2f\n', L/rd1)
fprintf('Near-far limit = %.2f\n', y_rel)
L_lmb = L/lambda;
limit_lmb = lim/lambda;

image_full1 = PSF_polar('N_rx', N_rx, 'y0', y0, 'bw',bw1,'c', c,'t_p', t_p,'d', d1,'fs',200e3, 'fc', fc, 'Nu', Nu); % [theta x range]
image_full2 = PSF_polar('N_rx', N_rx, 'y0', y0, 'bw',bw2,'c', c,'t_p', t_p,'d', d2,'fs',200e3, 'fc', fc, 'Nu', Nu); % [theta x range]

%%

psf1 = max(image_full1.data, [], 2);
psf2 = max(image_full2.data, [], 2);

test_area = trapezoid_area(N_rx, bw1, fc, c, image_full1.u, d1);
test_area2 = trapezoid_area(N_rx, bw2, fc, c, image_full2.u, d2);

f = figure('Position', [360 50.3333 769 590.6667]);

ax1=subplot(2,1,1);
hold on
plot(image_full1.u, db(abs(psf1./max(psf1))), 'LineWidth', 1);
plot(image_full2.u, db(abs(psf2./max(psf2))), 'LineWidth', 1);
set(ax1,'Xticklabel',[])
ylabel('Power [dB]')
xlim([image_full.u(1) image_full.u(end)])
ylim([-35 0])
title('Point spread function')
subtitle(sprintf('$%.f$ elements, $L/\\Delta r = %.1f$, $\\frac{y}{L_{max}^2/\\lambda} = %.2f$', N_rx, L/rd1, y_rel),'Interpreter','latex')
grid on

ax2=subplot(2,1,2);
hold on
plot(image_full1.u, test_area, 'LineWidth', 1)
plot(image_full2.u, test_area2,  'LineWidth', 1)
xlabel('$u$', 'Interpreter', 'latex')
ylabel('Area [m]')
title('Trapezoid area')

linkaxes([ax1,ax2],'x')
newcolors = [0/255 138/255 222/255;
            255/255 31/255 91/255;
             0/255 205/255 109/255;
             255/255 198/255 30/255;
             242/255 133/252 34/255];
colororder(newcolors)
grid on
set([ax1, ax2], 'LineWidth', 1, 'FontSize', 11, 'FontName', 'Serif')
legend(sprintf('$B_1/f_c = %.1f$,\n$d_1/\\lambda = %.1f$,\n$L_1 = %.f\\lambda$', bw1/fc, d1/lambda, L/lambda), ...
    sprintf('$B_2/f_c = %.1f$,\n$d_2/\\lambda = %.1f$,\n$L_2 = %.f\\lambda$', bw2/fc, d2/lambda, L2/lambda),'Interpreter', 'latex', 'Location', 'northeast')

% figname = sprintf('./Code/figs/grating-lobe-analysis/u-dependence/ex_rb1-%.2f_rb2-%g_d1-%g_d2-%g.png', bw1/fc, bw2/fc, d1/lambda, d2/lambda);
% exportgraphics(f, figname, 'Resolution', 300)
% figname = sprintf('./Code/figs/grating-lobe-analysis/u-dependence/ex_rb1-%.2f_rb2-%g_d1-%g_d2-%g.pdf', bw1/fc, bw2/fc, d1/lambda, d2/lambda);
% exportgraphics(f, figname, 'Resolution', 300)

