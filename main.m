clear all;
time = 60; %in seconds

folder = 'Data/denovo_2day/';
block_name.day2 = {'B1_baseline','B2_darkBaseline','B3','B5','B6_dark','B9','B10_dark','B11_habit'};
graph_name.day2 = {'Baseline','Baseline (D)','Early','Day 1','Day 1 (D)','Day 2','Day 2 (D)','Flip (F)'};
subj_name = {'subj1','subj3','subj4','subj5','subj6','subj7','subj8','subj9','subj10','subj11','subj12','subj13','subj14'};
remove = 2;
d = load_data(subj_name,block_name.day2,folder,time,remove);
data.day2 = analyze_data(d,block_name.day2);

folder = 'Data/denovo_5day/';
block_name.day5 = block_name.day2(1:end-3);
graph_name.day5 = graph_name.day2(1:end-1);
block_name.day5 = [block_name.day5 {'B10','B11_dark','B15','B16_dark','B20','B21_dark','B24','B25_dark','B26_habit'}];
graph_name.day5 = [graph_name.day5 {'Day 3','Day 3 (D)','Day 4','Day 4 (D)','Day 5','Day 5 (D)','Flip (F)'}];
subj_name = {'subj13','subj15','subj17','subj18','subj19','subj21','subj22','subj23','subj24','subj25','subj26','subj27','subj28','subj29'};
remove = 0;
d = load_data(subj_name,block_name.day5,folder,time,remove);
data.day5 = analyze_data(d,block_name.day5);


folder = 'Data/denovo_10day/';
block_name.day10 = block_name.day5(1:end-3);
graph_name.day10 = graph_name.day5(1:end-1);
block_name.day10 = [block_name.day10 {'B25','B26_dark','B30','B31_dark','B35','B36_dark','B40','B41_dark','B45','B46_dark','B49','B50_dark','B51_habit'}];
graph_name.day10 = [graph_name.day10 {'Day 6','Day 6 (D)','Day 7','Day 7 (D)','Day 8','Day 8 (D)','Day 9','Day 9 (D)','Day 10','Day 10 (D)','Flip (F)'}];
subj_name = {'subj1','subj2','subj3','subj4','subj5'};
remove = 1; % subj4 B1_baseline the Flock of Birds was bad
d = load_data(subj_name,block_name.day10,folder,time,remove);
data.day10 = analyze_data(d,block_name.day10);

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
% graph_coherence(data, 'cursor') % coherence
% graph_xcorr(data,block_name,graph_name)