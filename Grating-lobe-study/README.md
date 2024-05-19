## Dependencies

- Signal Processing Toolbox
- Symbolic Math Toolbox

## File overview

- **psf_gratinglobes_polar.m**: All the analyses of the grating lobes are done in this script. The script uses several functions that are described below. 
- **PSF_polar.m**: This function performs pulse-echo beamforming of a point scatterer. The function uses two functions to perform pulse-echo beamforming:
    - **pulsecomp.m**: Pulse compresses the received signal with the transmitted signal.
    - **BF_new.**: Performs beamforming. 
- **REL_LEVEL.m**: Calculates relative grating lobe level for all grating lobes in the point spread function.
- **trapezoid_area.m**: Calculates the trapezoid area of a grating lobe at u.
- **rel_blob_power.m**: Calculates the relative sum of power or the relative average power between mainlobe and grating lobe 1, grating lobe 2, etc. The function also has options to superpose trapezoid boxes on the point scatterer data and plot and save the images.