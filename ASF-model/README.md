# Broadband aperture smoothing function model

This folder provides the code for the broadband aperture smoothing function (ASF). 

## File overview

- **ASF_PW.m**: This function is the broadband ASF for a SAS with yaw error. It calculates the single frequency ASFs and sums all frequencies in the signal. The function can be run with one call without any inputs.
- **grating_lobe_BP.m**: The single-frequency ASF, used in the broadband ASF. It calculates the receiver beampattern, transmit beampattern, and element responses, and multiplies them to form the ASF product. The function can also plot the single-frequency ASF by providing the argument `plotting=1`. The weights on the receivers are uniform.
- **BP.m**: Calculates the modified transmit beampattern $W_{Tx,mod}(\theta)$.
