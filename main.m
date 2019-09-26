clear all;

time = 30; %in seconds
theta = 90;
groups = {'L','R'};
block_name = {'B1'};
subj_name = {'subj1','subj2','subj3','subj4'};
graph_name = {'baseline'};
folder = 'Data/stroke/left/';
d.L = load_data(subj_name,block_name,folder,time);
folder = 'Data/stroke/right/';
d.R = load_data(subj_name,block_name,folder,time);

data.L = analyze_data(d.L,block_name,false,0);
data.R = analyze_data(d.R,block_name,false,0);

% save dat data;
disp('Done')

%% graph phasor plots
subj = 'subj2';
graph_phasor(data, subj);

%% graph Bode plots: averaged across individuals
% average across individuals
% graph_bode(data, graph_name, groups);

%% simplified experimental data
gblocks = [2:5];
graph_bode_simple(data, graph_name, groups, gblocks);

% each subject individually
% subj_name2 = {'subj35'};
% graph_bode_solo(data,subj_name2,block_name,graph_name);

% for control2
% graph_bode2(data, graph_name, groups);

%% graph performance
gblocks = 1;

graph_amp_avg(data,groups,block_name,gblocks,graph_name); % amplitude spectrums
% graph_MSE(data, groups, block_name, graph_name); % mean squared error
% graph_rotatedMSE(data, groups, block_name, gblocks, graph_name);
% graph_angle(data, groups, gblocks,theta); % angle of movement error
% graph_lag(data, groups, gblocks, graph_name); % response lag
% graph_complexError(data, groups, graph_name, gblocks); % complex tracking error
% graph_coherence(data, groups, gblocks, graph_name) % coherence
