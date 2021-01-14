% clear all;
time = 40; %in seconds

folder = 'Data/online/pilot_learning';
remove = 0;
d = load_data(folder,time);
data = analyze_data(d);

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

graph_amp_avg(data,gblocks,subj,1); % amplitude spectrums
% graph_gainMatrix(data)
% graph_MSE(data, block_name, graph_name); % mean squared error
% graph_lag(data, gblocks, graph_name,'Lhand'); % response lag
% graph_complexError(data, graph_name, gblocks); % complex tracking error
% graph_coherence(data, block_name, graph_name, 'cursor') % coherence
% graph_xcorr(data,block_name,graph_name)