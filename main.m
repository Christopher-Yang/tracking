% VMR vs MR w/ no P2P training

clear all;
folder = 'Data/no_P2P/';
time = 40; %in seconds
theta = 90;
groups = {'rot','rot_i'};
block_name = {'baseline','rot1','rot2','rot3','rot4','after'};
graph_name = {'Baseline','Early','Train1','Train2','Late','Aftereffects'};
subj_name = {'subj1','subj2','subj3','subj4','subj5','subj6','subj7','subj8','subj9','subj10','subj11','subj12','subj13','subj14','subj15','subj16','subj17','subj18','subj19','subj20'};
subj_rot = {'subj2','subj5','subj7','subj8','subj9','subj10','subj13','subj15','subj17','subj19'};
subj_rot_i = {'subj1','subj3','subj4','subj6','subj11','subj12','subj14','subj16','subj18','subj20'};
d = load_data(subj_name,block_name,folder,time);
data.rot = analyze_data(d,subj_rot,block_name,false,0);
data.rot_i = analyze_data(d,subj_rot_i,block_name,false,0);

% save dat data;
disp('Done')

%% graph Bode plots: averaged across individuals
% average across individuals
% graph_bode(data, graph_name, groups);

%% simplified experimental data
gblocks = [1:2 5:6];
graph_bode_simple(data, graph_name, groups, gblocks);

%% graph performance
gblocks = [1:2 5:6];

graph_amp_avg(data,groups,block_name,gblocks,graph_name); % amplitude spectrums
% graph_MSE(data, groups, block_name, graph_name); % mean squared error
% graph_rotatedMSE(data, groups, block_name, gblocks, graph_name);
% graph_lag(data, groups, gblocks, graph_name); % response lag
% graph_complexError(data, groups, graph_name, gblocks); % complex tracking error
% graph_coherence(data, groups, gblocks, graph_name) % coherence
