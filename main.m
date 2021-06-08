clear all

% set variables for main experiment
folder = 'Data/fb_ff/'; % set path to data

% names of blocks: pert1=Early, pert4=Late, pert2 and pert3 are the
% tracking blacks that occurred between Early and Late
block_name = {'B1','B2','B3'}; 
subj_rot = {'subj1'}; % rotation subjects
subj_mir = {''}; % mirror-reversal subjects

% NOTE: subj1 had extra data at end of data files for each trial; need to
% go fix KineReach code

% Extract raw data. The first three fields in "d" are organized as 
% "d.[group]{[subject number]}.[block]". In the block subfield, "traj"
% contains the raw data for target, cursor, and hand position as well as
% timestamps of data collection. "tFile" contains the frequencies,
% amplitudes, and phases used to generate the target's sum-of-sinusoids
% trajectory.
d = load_data(subj_rot,subj_mir,block_name,folder);

% Perform primary data analysis on the rotation and mirror-reversal groups.
% analyze_data() extracts position and time information from the variable
% "d" and processes it to obtain Fourier transforms, mean-squared error, 
% and amplitude and phase information for the target, amongst other things.
% The first three fields of "data1" are organized as
% "data1.[group]{[subject number]}.[block]". For more information on the
% contents of the data structure, see analyze_data().
data.rot = analyze_data(d.rot,subj_rot,block_name);
% data.mir = analyze_data(d.mir,subj_mir,block_name);

% save data1 data1; % save data structure
disp('Done')

%% PART 2: FIGURES FOR MAIN EXPERIMENT
% load dat1 % load saved data

% Figure 2B: raw hand and target trajectories
graph_traj(data1);

% Figure 2C: mean-squared error
graph_MSE(data1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3: alignment matrix. "delay_opt" contains the optimal time to 
% align the target and hand trajectories that minimizes mean squared error.
% Loading delay_opt will significantly reduce computation time.
% If you would like to compute delay_opt instead of load it, comment 
% lines 72-73 and comment line 76. 

% USE PRECOMPUTED OPTIMAL DELAYS
load delay_opt
graph_alignMatrix(data1,delay_opt);

% COMPUTE OPTIMAL DELAYS FROM SCRATCH (takes ~60 mins to run)
% graph_alignMatrix(data1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure 4A and Figure 4-supplements 1A & 2: amplitude spectra
graph_ampSpectra(data1);

% Figures 4B and Figure 4-supplement 1B: coherence
graph_coherence(data1);

% Figures 4C and Figure 4-supplement 1C: lag 
graph_lag(data1);

% Figure 5, Figure 5-supplements 1 & 2: gain matrix
graph_gainMatrix(data1,experiment);
