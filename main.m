% clear all;
time = 40; %in seconds

folder = 'Data/online/pilot_learning2';
subj_rot = {'pilot13','pilot14','pilot15','pilot16','pilot17'};
subj_mir = {'pilot8','pilot9','pilot10','pilot12'};
remove = 0;
d.rot = load_data(folder,subj_rot,time);
d.mir = load_data(folder,subj_mir,time);

data.rot = analyze_data(d.rot);
data.mir = analyze_data(d.mir);

% save dat data;
disp('Done')

%% graph Bode plots: averaged across individuals
% average across individuals
% graph_bode(data, graph_name);

%% simplified experimental data
gblocks = 1:3;
graph_bode_simple(data, graph_name, gblocks,'cursor');

%% graph performance
gblocks = 1;
subj = [];

graph_amp_avg(data,gblocks,1,1); % amplitude spectrums
% graph_gainMatrix(data)
% graph_MSE(data, block_name, graph_name); % mean squared error
% graph_lag(data, gblocks, graph_name,'Lhand'); % response lag
% graph_complexError(data, graph_name, gblocks); % complex tracking error
% graph_coherence(data, block_name, graph_name, 'cursor') % coherence
% graph_xcorr(data,block_name,graph_name)