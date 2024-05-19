function grels = REL_LEVEL(psf, lambda, d, c, fc, bw, u, L, N_rx)
% Calculates relative grating lobe level for all grating lobes
% assumes grating lobes are symmetric.
% IN:
%       psf (complex array) : Point spread function
% OUT:
%       grels : Relative grating lobe levels, first element is p=1.

nr_globes = round(d/lambda) - 1; % Nr. of grating lobes in visible region
if d/lambda > 1 && d/lambda < 1.5
    nr_globes = ceil(d/lambda) - 1;
end

grels = zeros(nr_globes, 1); 

% [mainlobe_val,loc] = max(abs(psf)); % Main lobe index
loc = round(length(psf)/2);
mainlobe_val = abs(psf(loc));
bw_globew = c/d*(1/(fc-bw/2)-1/(fc+bw/2));
du = u(2)-u(1);

for k=1:length(grels)
%     globe_val = abs(psf(loc + round(k*lambda/d/du))); % Just calculate
%     the max. rel level at fc 
%     grels(k) = db(globe_val/mainlobe_val);
% end
    if (2*lambda/L) < bw_globew
        % If the width due to broadband is larger than mainlobe width 
        glint = (loc+round(k*c/(fc+bw/2)/d/du)):(loc+round(k*c/(fc-bw/2)/d/du)); % Interval for grating lobe i

        if (glint(end) > length(psf))
            glint = (loc+round(k*c/(fc+bw/2)/d/du)):length(psf);
        end
        globe_val = max(abs(psf(glint)));

        grels(k) = db(globe_val/mainlobe_val);
    else
        glint = (loc + round((k-1/N_rx)*lambda/(d*du))):(loc + round((k+1/N_rx)*lambda/(d*du)));

        if (glint(end) > length(psf))
            glint = (loc + round((k-1/N_rx)*lambda/(d*du))):length(psf);
        end
        globe_val = max(abs(psf(glint)));

        grels(k) = db(globe_val/mainlobe_val); % Interval for grating lobe i

    end
end