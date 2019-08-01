clear all;

folder = 'Data/bimanual_pilot/';
time = 40; %in seconds
theta = 90;
groups = {'bim'};
% block_name = {'B1','B2','B8','B21','B22_habit'};
% graph_name = {'Baseline','Early','Middle','Late','Habit'};
% subj_name = {'subj1','subj2','subj3'};
block_name = {'D1/B1_baseline','D1/B4','D2/B10','D3/B16','D4/B22','D5/B28'};
graph_name = {'B1','B4','B10','B16','B22','B28'};
subj_name = {'subj4'};
d = load_data(subj_name,block_name,folder,time);
block_name = regexprep(block_name,'^...','');
data = analyze_data(d,block_name,0);

% save dat data;
disp('Done')

%% graph phasor plots
graph_name2 = {'Baseline', 'Naive', 'Train1','Train2','Max training', 'Aftereffects'};
graph_phasor(data, block_name, groups, graph_name2, subj_rot, subj_rot_i);

%% graph Bode plots: averaged across individuals
% average across individuals
% graph_bode(data, graph_name, groups);

%% simplified experimental data
gblocks = [2:5];
graph_bode_simple(data, graph_name, groups, gblocks);

%% graph performance
gblocks = 1:2;

graph_amp_avg(data,groups,block_name,gblocks,graph_name); % amplitude spectrums
% graph_MSE(data, groups, block_name, graph_name); % mean squared error
% graph_rotatedMSE(data, groups, block_name, gblocks, graph_name);
% graph_angle(data, groups, gblocks,theta); % angle of movement error
% graph_lag(data, groups, gblocks, graph_name); % response lag
% graph_complexError(data, groups, graph_name, gblocks); % complex tracking error
% graph_coherence(data, groups, gblocks, graph_name) % coherence
