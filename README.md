# Figures for tracking task #

This collection is associated with the article "De novo learning
versus adaptation of continuous control in a manual tracking task."
This code generates all figures related to the sum-of-sinusoids
tracking task. Raw data is contained in the folder "Data." Data from
the main experiment (used to generate Figures 2B-5C and associated 
figure supplements) is contained in "Data/vmr90_vs_mr" while data from 
the second experiment (used to generate Figure 6 and associated figure 
supplements) is contained in "Data/no_p2p." Each experiment's data is 
first organized by subject, then block. The block name corresponds to 
the following:

    baseline: tracking without the rotation or mirror-reversal
    pert1: first block of tracking under either perturbation
    pert2: second block "..."
    pert3: third block "..."
    pert4: fourth block "..."
    post: tracking without the perturbations to assess aftereffects

Data from each trial is stored as a separate .dat file within each
block folder, labeled as "Data*.dat". The columns in each .dat file
correspond to the following:

    columns 1-2: target x/y position
    columns 3-4: left hand x/y position (not used for this paper)
    columns 5-6: right hand x/y position
    columns 7-8: cursor x/y position
    columns 9-10: information not used for this paper
    column 11: timestamp of data collection

Also stored is a "tFile" which defines the frequencies,
amplitudes, and phases used to generate the target's sum-of-sinusoids
trajectory. The values in the tFile can be read as follows:

    values 1-7: x frequencies
    values 8-14: y frequencies
    values 15-21: x amplitudes
    values 22-28: y amplitudes
    values 29-35: x phases
    values 36-42: y phases

All analyses can be performed and all figures generated by simply
running main.m. All other .m files are functions that are used by
main.m for analysis or plotting. Briefly, these functions do the
following:

    analyze_data.m: primary data analysis
    editErrorBar.m: edits plots generated by shadedErrorBar.m
    error_ellipse.m: plots confidence ellipse
    fourier.m: discrete Fourier transform and phasor analysis
    graph_alignMatrix.m: plots alignment matrices
    graph_ampSpectra.m: plots amplitude spectra
    graph_coherence.m: plots spectral coherence plots
    graph_gainMatrix.m: plots gain matrices
    graph_MSE.m: plots mean-squared error
    graph_traj.m: plots target and cursor trajectories
    load_data.m: extracts raw data from data files
    LQR.m: simulates a time-based catch-up strategy model
    LQR_distance.m: a distance-based catch-up strategy model
    shadedErrorBar.m: plots lines with shaded error bars

The following *.mat files are included that contain precomputed
matrices that dramatically reduce computation time (see code for
details):

    delay_opt.mat: used by graph_alignMatrix.m
    simResults.mat: used by LQR.m

This code requires the following MATLAB Toolboxes:

    Optimization Toolbox
    DSP System Toolbox
    Signal Processing Toolbox
    Phased Array System Toolbox

Two scripts included here were created by other individuals,
shadedErrorBar.m and error_ellipse.m. The associated licenses for each
script can be found in the Licenses folder.