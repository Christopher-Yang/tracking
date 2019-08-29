clear all;

folder = 'Data/bimanual_pilot/';
time = 60; %in seconds
theta = 90;
% block_name = {'B1','B2','B8','B21','B22_habit'};
% graph_name = {'Baseline','Early','Middle','Late','Habit'};
% subj_name = {'subj1','subj2','subj3'};

% block_name = {'B1_baseline','B11','B14','B17','B20','B23'};
% graph_name = {'B1','B11','B14','B17','B20','B23'};
% subj_name = {'subj5','subj6'};

block_name = {'B1_baseline','B5','B16','B26'};
graph_name = {'B1','B5','B16','B26'};
subj_name = {'subj7','subj8','subj9'};
d = load_data(subj_name,block_name,folder,time);
data = analyze_data(d,block_name,0);

% save dat data;
disp('Done')

%% graph Bode plots: averaged across individuals
% average across individuals
% graph_bode(data, graph_name);

%% simplified experimental data
gblocks = 1:4;
graph_bode_simple(data, graph_name, gblocks);

%% graph performance
gblocks = 1:6;

% graph_amp_avg(data,block_name,gblocks,graph_name); % amplitude spectrums
graph_MSE(data, block_name, graph_name); % mean squared error
% graph_rotatedMSE(data, block_name, gblocks, graph_name);
% graph_lag(data, gblocks, graph_name); % response lag
% graph_complexError(data, graph_name, gblocks); % complex tracking error
% graph_coherence(data, gblocks, graph_name) % coherence
