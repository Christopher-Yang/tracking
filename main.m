clear all;

% subj 1: me
% subj 2: Martina
% subj 3: me
% subj 4: Kahori
% subj 5: Yue
% subj 6: Alex
% subj 7: Tim
% subj 8 onward: other participants

folder = 'Data/baseline_learning/';
time = 60; %in seconds
block_name = {'B1','B2','B3','B4','B5','B6','B7'};
graph_name = block_name;
% block_name = {'right1','F1_1','F3_1','F4_1'};
% graph_name = block_name;
subj_name = {'subj8','subj9','subj10','subj11'};
d = load_data(subj_name,block_name,folder,time);
data = analyze_data(d,block_name,0,0);

% save dat data;
disp('Done')

%% graph Bode plots: averaged across individuals
% average across individuals
% graph_bode(data, graph_name);

%% simplified experimental data
gblocks = [2:5];
graph_bode_simple(data, graph_name, gblocks);

%% graph performance
gblocks = 1:7;
output = 'Rhand';

% graph_amp_avg(data,block_name,gblocks,graph_name,output); % amplitude spectrums
% graph_MSE(data, block_name, graph_name); % mean squared error
% graph_rotatedMSE(data, block_name, gblocks, graph_name);
% graph_angle(data, gblocks,theta); % angle of movement error
graph_lag(data, gblocks, graph_name,'cursor'); % response lag
% graph_complexError(data, graph_name, gblocks); % complex tracking error
% graph_coherence(data, gblocks, graph_name) % coherence
