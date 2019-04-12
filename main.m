clear all;

folder = 'Data/bimanual_pilot/';
time = 40; %in seconds
theta = 90;
groups = {'bim'};
block_name = {'B1','B2','B8','B21'};
graph_name = {'Baseline','Early','Middle','Late'};
subj_name = {'subj1','subj2','subj3'};
d = load_data(subj_name,block_name,folder,time);
data.bim = analyze_data(d,subj_name,block_name,false,0);

% save dat data;
disp('Done')

%% graph phasor plots
graph_name2 = {'Baseline', 'Naive', 'Train1','Train2','Max training', 'Aftereffects'};
graph_phasor(data, block_name, groups, graph_name2, subj_rot, subj_rot_i);

%% graph Bode plots: averaged across individuals
% average across individuals
% graph_bode(data, graph_name, groups);

%% simplified experimental data
gblocks = [1:4];
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
gblocks = [1:4];

graph_amp_avg(data,groups,block_name,gblocks,graph_name); % amplitude spectrums
% graph_MSE(data, groups, block_name, graph_name); % mean squared error
% graph_rotatedMSE(data, groups, block_name, gblocks, graph_name);
% graph_angle(data, groups, gblocks,theta); % angle of movement error
% graph_lag(data, groups, gblocks, graph_name); % response lag
% graph_complexError(data, groups, graph_name, gblocks); % complex tracking error
% graph_coherence(data, groups, gblocks, graph_name) % coherence
