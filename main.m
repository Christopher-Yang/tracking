% VMR90 with no cursor feedback during tracking but with feedback during 150 P2P

clear all;

folder = 'Data/dark/';
time = 40; %in seconds
theta = 90;
groups = {'rot'};
block_name = {'baseline','dark_baseline','rot1','rot2','rot3','rot4','after'};
graph_name = {'Baseline','Dark Baseline','Early','rot2','rot3','Late','Post'};
subj_name = {'subj1','subj2','subj3','subj4','subj5'};
subj_rot = subj_name;
d = load_data(subj_name,block_name,folder,time);
data = analyze_data(d,block_name,0);

% save dat data;
disp('Done')

%% graph Bode plots: averaged across individuals
% average across individuals
% graph_bode(data, graph_name, groups);

%% simplified experimental data
gblocks = [1:2 5:6];
graph_bode_simple(data, graph_name, groups, gblocks);

%% graph performance
gblocks = [1:2 6:7];

graph_amp_avg(data,block_name,gblocks,graph_name,'Rhand'); % amplitude spectrums
% graph_MSE(data, block_name, graph_name); % mean squared error
% graph_rotatedMSE(data, block_name, gblocks, graph_name);
% graph_lag(data, gblocks, graph_name, 'Rhand'); % response lag
% graph_complexError(data, graph_name, gblocks); % complex tracking error
% graph_coherence(data, block_name, gblocks, graph_name, 'Rhand') % coherence
