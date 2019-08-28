clear all;

folder = 'Data/baseline_learning/';
time = 60; %in seconds
block_name = {'left1','left4','right1','right4'};
graph_name = {'Left Early','Left Late','Right Early','Right Late'};
subj_name = {'subj1'};
d = load_data(subj_name,block_name,folder,time);
data = analyze_data(d,block_name,0);

% save dat data;
disp('Done')

%% graph Bode plots: averaged across individuals
% average across individuals
% graph_bode(data, graph_name);

%% simplified experimental data
gblocks = [2:5];
graph_bode_simple(data, graph_name, gblocks);

%% graph performance
gblocks = 1:4;

% graph_amp_avg(data,block_name,gblocks,graph_name); % amplitude spectrums
% graph_MSE(data, block_name, graph_name); % mean squared error
% graph_rotatedMSE(data, block_name, gblocks, graph_name);
% graph_angle(data, gblocks,theta); % angle of movement error
graph_lag(data, gblocks, graph_name,'cursor'); % response lag
% graph_complexError(data, graph_name, gblocks); % complex tracking error
% graph_coherence(data, gblocks, graph_name) % coherence
