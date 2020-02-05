% This file is associated with the article "De novo learning versus
% adaptation of continuous control in a manual tracking task." Running
% main.m generates all figures related to the sum-of-sinusoids tracking 
% task. 
% 
% Parts 1-2 of the script (lines 19-65) analyzes and makes figures for the 
% data from the main experiment. Parts 3-4 (lines 66-87) analyzes and 
% makes figures for the second experiment. Analyses and figures can be run
% independently of one another.
% 
% To generate all figures, run the entire script, which should take 
% between 5-10 minutes. If only certain figures are needed, first run the 
% analysis in Part 1 or Part 3. Then run the appropriate functions to 
% generate figures in Part 2 or Part 4, respectively. To avoid running the
% analyses multiple times, the variables "data1" and "data2", 
% corresponding to the data from the main and second experiment 
% respectively, can be saved and loaded to generate figures as needed.

%% PART 1: ANALYSIS FOR MAIN EXPERIMENT
clear all;

% set variables for analysis
experiment = 1; % 1: main experiment with point-to-point practice; 2: second experiment with minimal point-to-point
folder = 'Data/VMR90_vs_MR/'; % set path to data
block_name = {'baseline','pert1','pert2','pert3','pert4','post'}; % names of directories containing data
subj_rot = {'subj17','subj18','subj21','subj22','subj24','subj25','subj28','subj31','subj32','subj33'}; % rotation subjects
subj_mir = {'subj14','subj15','subj16','subj19','subj23','subj26','subj27','subj29','subj30','subj34'}; % mirror-reversal subjects

% analysis
d = load_data(subj_rot,subj_mir,block_name,folder); % extract raw data and store in d
data1.rot = analyze_data(d.rot,subj_rot,block_name,1); % analyze rotation group's data and store output in the variable data
data1.mir = analyze_data(d.mir,subj_mir,block_name,2); % analyze mirror-reversal group's data

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
graph_amp_avg(data1);

% Figures 4B and S1B: coherence
graph_coherence(data1);

% Figures 5, S2B, and S4: gain matrix
graph_gainMatrix(data1,experiment); 

% Figure S3: phasors
% rotation group
graph_phasor(data1.rot{4}.baseline.phasors.Rhand,22) % baseline
graph_phasor(data1.rot{4}.pert4.phasors.Rhand,23) % late learning
% mirror-reversal group
graph_phasor(data1.mir{9}.baseline.phasors.Rhand,24) % baseline
graph_phasor(data1.mir{9}.pert4.phasors.Rhand,25) % late learning

%% PART 3: ANALYSIS FOR SECOND EXPERIMENT
% set variables for analysis
experiment = 2;
folder = 'Data/no_P2P/'; % set path to data
time = 40; % seconds of data to be analyzed per trial; can be 20 or 40 secs
block_name = {'baseline','pert1','pert2','pert3','pert4','post'}; % names of directories containing data
subj_rot = {'subj2','subj5','subj7','subj8','subj9','subj10','subj13','subj15','subj17','subj19'}; % rotation subjects
subj_mir = {'subj1','subj3','subj4','subj6','subj11','subj12','subj14','subj16','subj18','subj20'}; % mirror-reversal subjects

% analysis
d = load_data(subj_rot,subj_mir,block_name,folder); % extract raw data and store in d
data2.rot = analyze_data(d.rot,subj_rot,block_name,1); % analyze rotation group's data and store output in the variable data
data2.mir = analyze_data(d.mir,subj_mir,block_name,2); % analyze mirror-reversal group's data

% save dat2 data2; % save data structure
disp('Done')

%% PART 4: FIGURES FOR SECOND EXPERIMENT
% load dat2 % load saved data

% Figure S5: gain matrix
graph_gainMatrix(data2,experiment); 