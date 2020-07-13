% clear all;
time = 40; %in seconds

folder = 'Data/online/';
% block_name = {'baseline','cursor_sines'};
% graph_name = {'Baseline','Cursor Sines'};
subj_name = {'Chris'};
remove = 0;
d = load_data(subj_name,folder,time,remove);
% data = analyze_data(d);

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

graph_amp_avg(data,block_name,gblocks,graph_name,'cursor'); % amplitude spectrums
% graph_MSE(data, block_name, graph_name); % mean squared error
% graph_lag(data, gblocks, graph_name,'Lhand'); % response lag
% graph_complexError(data, graph_name, gblocks); % complex tracking error
% graph_coherence(data, block_name, graph_name, 'cursor') % coherence
% graph_xcorr(data,block_name,graph_name)