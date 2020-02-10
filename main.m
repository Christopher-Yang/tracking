% This file is associated with the article "De novo learning versus
% adaptation of continuous control in a manual tracking task." Running
% main.m generates all figures related to the sum-of-sinusoids tracking 
% task. 
% 
% Parts 1-2 of the script (lines 19-83) analyzes and makes figures for the 
% data from the main experiment (Figures 2B-5 and S1-S4). Parts 3-4 (lines 
% 84-111) analyzes and makes figures for the second experiment (Figure S5). 
% Analysis and plotting code can be run independently of one another.
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

% save dat1 data1; % save data structure
disp('Done')

%% PART 2: FIGURES FOR MAIN EXPERIMENT
% load dat1 % load saved data

% Figure 2B: raw hand and target trajectories
graph_traj(data1); 

% Figure 2C: mean-squared error
graph_MSE(data1);

% Figure 3: transformation matrix
graph_transformMat(data1);

% Figures 4A, S1A, and S2A: amplitude spectra
graph_ampSpectra(data1);

% Figures 4B and S1B: coherence
graph_coherence(data1);

% Figures 5, S2B, and S4: gain matrix
graph_gainMatrix(data1,experiment); 

% Figure S3: phasors
% rotation group
graph_phasor(data1.rot{4}.baseline.Rhand.phasors,22) % baseline
graph_phasor(data1.rot{4}.pert4.Rhand.phasors,23) % late learning
% mirror-reversal group
graph_phasor(data1.mir{9}.baseline.Rhand.phasors,24) % baseline
graph_phasor(data1.mir{9}.pert4.Rhand.phasors,25) % late learning

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

% save dat2 data2; % save data structure
disp('Done')

%% PART 4: FIGURES FOR SECOND EXPERIMENT
% load dat2 % load saved data

% Figure S5: gain matrix
graph_gainMatrix(data2,experiment); 