% This file is associated with the article "De novo learning versus
% adaptation of continuous control in a manual tracking task." Running
% main.m generates all figures related to the sum-of-sinusoids tracking 
% task. 
% 
% Parts 1-2 of the script analyzes and makes figures for the data from the 
% main experiment (Figures 2B-5C and associated figure supplements). Parts 
% 3-4 analyzes and makes figures for the second experiment (Figure 6 and 
% associated figure supplements).
%
% To generate all figures, run the entire script, which should take 
% between 5-10 minutes. If only certain figures are needed, first run the 
% analysis in Part 1 or Part 3. Then run the appropriate functions to 
% generate figures in Part 2 or Part 4, respectively. To avoid running the
% analyses multiple times, the variables "data1" and "data2", 
% corresponding to the data from the main and second experiment 
% respectively, can be saved and loaded to generate figures as needed.

%% PART 1: ANALYSIS FOR MAIN EXPERIMENT
clear all

% set variables for main experiment
experiment = 1; % main experiment
folder = 'Data/vmr90_vs_mr/'; % set path to data

% names of blocks: pert1=Early, pert4=Late, pert2 and pert3 are the
% tracking blacks that occurred between Early and Late
block_name = {'baseline','pert1','pert2','pert3','pert4','post'}; 
subj_rot = {'subj17','subj18','subj21','subj22','subj24','subj25',...
    'subj28','subj31','subj32','subj33'}; % rotation subjects
subj_mir = {'subj14','subj15','subj16','subj19','subj23','subj26',...
    'subj27','subj29','subj30','subj34'}; % mirror-reversal subjects

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
data1.rot = analyze_data(d.rot,subj_rot,block_name);
data1.mir = analyze_data(d.mir,subj_mir,block_name);

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

%% PART 3: ANALYSIS FOR SECOND EXPERIMENT
% set variables for second experiment
experiment = 2; % second experiment
folder = 'Data/no_p2p/'; % set path to data

% names of blocks: pert1=Early, pert4=Late, pert2 and pert3 are the
% tracking blacks that occurred between Early and Late
block_name = {'baseline','pert1','pert2','pert3','pert4','post'};
subj_rot = {'subj2','subj5','subj7','subj8','subj9','subj10','subj13',...
    'subj15','subj17','subj19'}; % rotation subjects
subj_mir = {'subj1','subj3','subj4','subj6','subj11','subj12','subj14',...
    'subj16','subj18','subj20'}; % mirror-reversal subjects

% extract raw data; this performs the same function as load_data in Part 1
d = load_data(subj_rot,subj_mir,block_name,folder);

% analyze data; this performs the same function as analyze_data in Part 1
data2.rot = analyze_data(d.rot,subj_rot,block_name);
data2.mir = analyze_data(d.mir,subj_mir,block_name);

% save data2 data2; % save data structure
disp('Done')

%% PART 4: FIGURES FOR SECOND EXPERIMENT
% load dat2 % load saved data

% Figure 6, and Figure 6-supplement 1: gain matrix
graph_gainMatrix(data2,experiment);