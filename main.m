clear all;

folder = 'Data/bimanual_pilot/';
time = 40; %in seconds

% block_name = {'B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11','B12','B13','B14','B15','B16','B17','B18','B19','B21'};
% graph_name = block_name;
% subj_name = {'subj1','subj2','subj3'};
% d = load_data(subj_name,block_name,folder,time);
% rotate = 0;
% data = analyze_data(d,block_name,0,rotate);

% subj 4 had trial lengths too short; subj 5 was poor at task
% block_name = {'B1_baseline','B4','B6','B8','B11','B14','B17','B20','B23'};
% graph_name = {'Baseline','B4','B6','B8','B11','B14','B17','B20','B23'};
% subj_name = {'subj6'};
% d = load_data(subj_name,block_name,folder,time);
% rotate = 0;
% data = analyze_data(d,block_name,0,rotate);

% block_name = {'B1_baseline','B4_darkBaseline','B28_dark'};
% graph_name = {'Baseline','Dark baseline','Dark Late'};
% block_name = {'B1_baseline','B5','B6','B7','B8','B9','B10','B13','B14','B15','B16','B19','B20','B21','B22','B25','B26'};
% graph_name = {'Baseline','Early','Late'};
% subj_name = {'subj7','subj8','subj9'};
% d = load_data(subj_name,block_name,folder,time);
% rotate = 0;
% data = analyze_data(d,block_name,0,rotate);

block_name = {'B1_baseline','B5','B6','B9','B10','B11','B12','B15','B16','B17','B18','B22','B23','B24','B27','B28','B29'};
graph_name = {'B1','B5','B6','B9','B10','B11','B12','B15','B16','B17','B18','B22','B23','B24','B27','B28','B29'};
% block_name = {'B1_baseline','B3_dualBaseline','B5','B6','B7_dual','B9','B10','B11','B12','B13_dual','B15','B16','B17','B18','B19_dual','B22','B23','B24','B25_dual','B27','B28','B29','B30_dual','B32_habit','B33_habit'};
% graph_name = {'B1','B3 (d)','B5','B6','B7 (d)','B9','B10','B11','B12','B13 (d)','B15','B16','B17','B18','B19 (d)','B22','B23','B24','B25 (d)','B27','B28','B29','B30 (d)','B32 (h)','B33 (h)'};
subj_name = {'subj10','subj11'};
d = load_data(subj_name,block_name,folder,time);
rotate = 1;
data = analyze_data(d,block_name,0,rotate);

% save dat data;
disp('Done')

%% graph Bode plots: averaged across individuals
% average across individuals
% graph_bode(data, graph_name);

%% simplified experimental data
gblocks = 1:3;
graph_bode_simple(data, graph_name, gblocks,'cursor');

%% graph performance
gblocks = [1 4 10 20];
subj = 1:2;

% graph_amp_avg(data,block_name,gblocks,graph_name,'Rhand'); % amplitude spectrums
graph_MSE(data, block_name, graph_name); % mean squared error
% graph_lag(data, gblocks, graph_name,'cursor'); % response lag
% graph_complexError(data, graph_name, gblocks); % complex tracking error
% graph_coherence(data, block_name, graph_name, 'cursor', subj) % coherence
% graph_xcorr(data,block_name,graph_name)