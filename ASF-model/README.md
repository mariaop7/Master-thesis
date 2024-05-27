# Broadband aperture smoothing function model

This folder provides the code for the broadband aperture smoothing function (ASF) and an example of how to run it.

## File overview

- **ASF_PW.m**: This function is the broadband ASF for a SAS with yaw error. It calculates the single frequency ASFs and sums all frequencies in the signal. The function can be run with one call without any inputs as seen in the example below.
    - **grating_lobe_BP.m**: The single-frequency ASF used in the broadband ASF. It calculates the receiver beampattern, transmit beampattern, and element responses, and multiplies them to form the ASF product. The function can also plot the single-frequency ASF by providing the argument `plotting=1`. The weights on the receivers are assumed to be uniform.
    - **BP.m**: Calculates the modified transmit beampattern $W_{Tx,mod}(\theta)$.

<details>
<summary>Example</summary>

```Matlab
[total_BP, angles, u] = ASF_PW(); % Create broadband ASF with pre-defined parameters

f = figure('Position', [314 135 766 362]);

plot(u, db(abs(total_BP/max(total_BP(:)))), 'LineWidth', 1)
grid on
xlim([-0.1 0.1])
ylim([-75 5])
xlabel('$u$', 'Interpreter','latex')
ylabel('Power [dB]')
title('Periodic error grating lobes')
subtitle('Aperture smoothing function model')
set(gca,'LineWidth', 1, 'Fontsize', 12, 'FontName', 'Serif')
```
</details>
