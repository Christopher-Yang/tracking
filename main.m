clear all;

folder = 'Data/VMR90_vs_MR/';
time = 40; %in seconds
theta = 90;
groups = {'rot','rot_i'};
block_name = {'no_rot1','rot1','rot2','rot3','rot4','no_rot2'};
graph_name = {'Baseline','Early','Train2','Train3','Late','Post'};
subj_name = {'subj14','subj15','subj16','subj17','subj18','subj19','subj21','subj22','subj23','subj24','subj25','subj26','subj27','subj28','subj29','subj30','subj31','subj32','subj33','subj34'};
subj_rot = {'subj17','subj18','subj21','subj22','subj24','subj25','subj28','subj31','subj32','subj33'};
subj_rot_i = {'subj14','subj15','subj16','subj19','subj23','subj26','subj27','subj29','subj30','subj34'};
d = load_data(subj_name,block_name,folder,time);
data.rot = analyze_data(d,subj_rot,block_name,false,1);
data.rot_i = analyze_data(d,subj_rot_i,block_name,false,2);

% save dat data;
disp('Done')

%% graph phasor plots
graph_name2 = {'Baseline', 'Naive', 'Train1','Train2','Max training', 'Aftereffects'};
graph_phasor(data, block_name, groups, graph_name2, subj_rot, subj_rot_i);

%% graph Bode plots: averaged across individuals
% average across individuals
% graph_bode(data, graph_name, groups);

%% simplified experimental data
gblocks = [1:2 5:6];
graph_bode_simple(data, graph_name, groups, gblocks);
% graph_bode_hilo(data, groups, graph_name,1:6);

% simplified control1 
% graph_name2 = {'Baseline', 'No training', 'Max training', 'Aftereffects', 'No training (opp)', 'Opp train 1'};
% gblocks = [1:2 5:7];
% graph_bode_hilo2(data,groups,graph_name,gblocks);

% each subject individually
% subj_name2 = {'subj35'};
% graph_bode_solo(data,subj_name2,block_name,graph_name);

% for control2
% graph_bode2(data, graph_name, groups);

%% graph performance
gblocks = [1:2 5:6];

graph_amp_avg(data,groups,block_name,gblocks,graph_name); % amplitude spectrums
% graph_MSE(data, groups, block_name, graph_name); % mean squared error
% graph_rotatedMSE(data, groups, block_name, gblocks, graph_name); % rotating MSE
% graph_angle(data, groups, gblocks,theta); % angle of movement error
% graph_lag(data, groups, gblocks, graph_name); % response lag
% graph_complexError(data, groups, graph_name, gblocks); % complex tracking error
% graph_coherence(data, groups, gblocks, graph_name) % coherence
