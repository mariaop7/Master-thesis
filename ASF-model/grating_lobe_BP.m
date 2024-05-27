function [W_2w, kx] = grating_lobe_BP(N, M, skew_angle, L_tx, d_Rx, d_Tx, c, Fc, u, plotting, w_Tx)  
% Returns Tx and Rx beampattern of the SAS at single frequency Fc. Assumes
% uniform receiver weights.
% IN:
%       N               Receive elements
%       M               Transmits
%       skew_angle      Yaw error
%       L_tx            Distance between pings 
%       d_Rx            Interelement distance in receiver array
%       d_Tx            Transmitter element width
%       c               Speed of sound
%       Fc              Frequency
%       u               u-space values of ASF
%       plotting        Whether produce plots (plotting=1)
%       w_Tx            Tapering on the transmitters
% OUT:
%       W_2w            The total ASF
%       kx              Corresponding wavenumber kx values

lambda = c/Fc; % Wavelength

k0 = 2*pi/lambda*sin(skew_angle/180*pi);
kx = 2*pi/lambda*u; % Wavenumber

% Receive beampattern
W_sub_Rx = 1/N*sin((kx+k0)*d_Rx*N/2)./sin((kx+k0)*d_Rx/2); % For uniform receiver weights

% Modified transmit array + beampattern
ELpos_Tx = (0:(M-1)).*L_tx - (M-1)*L_tx/2; % Assumes centered transmit array
W_Tx_mod = BP(ELpos_Tx,kx*2, w_Tx); 

% Element responses
W_el_Tx = sin(kx*d_Tx/2)./(kx*d_Tx/2); 
W_el_Rx = sin(kx*d_Rx/2)./(kx*d_Rx/2);

% Total ASF
W_2w = W_Tx_mod.*W_sub_Rx.*W_el_Tx.*W_el_Rx;

% Whether plotting
if plotting
    plot(u, db(abs(W_Tx_mod)),'LineWidth', 1);
    hold on 
    plot(u, db(abs(W_sub_Rx)), 'LineWidth', 1);
    plot(u, db(abs(W_2w)), 'LineWidth',  1);
    hold off
    xlabel('$u$', 'Interpreter', 'latex')
    ylabel('Power [dB]')
    ylim([-50 5])
    title('Aperture smoothing function')
    subtitle(sprintf('Yaw error $= %.2f^\\circ$, $Tx = %2d$, $Rx = %2d$, $L_{Tx} = %g$ m', skew_angle, M, N, L_tx), 'Interpreter','latex')
    legend('W_{Tx,mod}','W_{Rx}(\theta+\theta_0)','W_{tot}')
    set(gca, 'LineWidth', 1)
    grid on 

    newcolors = [0/255 205/255 109/255;
                 0/255 138/255 222/255;
                 255/255 31/255 91/255;
                 255/255 198/255 30/255;
                 242/255 133/252 34/255];
    colororder(newcolors)
end