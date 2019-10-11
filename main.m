clear all;

folder = 'Data/VMR90_vs_MR/';
time = 40; %in seconds
theta = 90;
groups = {'rot','mir'};
block_name = {'no_rot1','rot1','rot2','rot3','rot4','no_rot2'};
graph_name = {'Baseline','Early','Train2','Train3','Late','Post'};
subj_rot = {'subj17','subj18','subj21','subj22','subj24','subj25','subj28','subj31','subj32','subj33'};
subj_mir = {'subj14','subj15','subj16','subj19','subj23','subj26','subj27','subj29','subj30','subj34'};
d = load_data(subj_rot,subj_mir,block_name,folder,time);
data.rot = analyze_data(d.rot,subj_rot,block_name,1);
data.mir = analyze_data(d.mir,subj_mir,block_name,2);

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

graph_amp_avg(data,groups,block_name,gblocks,graph_name,'Rhand'); % amplitude spectrums
% graph_MSE(data, groups, block_name, graph_name); % mean squared error
% graph_rotatedMSE(data, groups, block_name, gblocks, graph_name); % rotating MSE
% graph_lag(data, groups, gblocks, graph_name); % response lag
% graph_complexError(data, groups, graph_name, gblocks); % complex tracking error
% graph_coherence(data, groups, block_name, gblocks, graph_name,'Rhand') % coherence
