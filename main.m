clear all;

folder = 'Data/bimanual_pilot/';

% time = 40; %in seconds
% block_name = {'B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11','B12','B13','B14','B15','B16','B17','B18','B19','B21','B22_habit','B23','B24'};
% graph_name = {'B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11','B12','B13','B14','B15','B16','B17','B18','B19','B21','B22(H)','B23','B24'};
% subj_name = {'subj1','subj2','subj3'};
% d = load_data(subj_name,block_name,folder,time);
% rotate = 0;
% data = analyze_data(d,block_name,0,rotate);

% subj 4 had trial lengths too short; subj 5 was poor at task
% time = 60; %in seconds
% block_name = {'B1_baseline','B4','B6','B8','B11','B14','B17','B20','B23'};
% graph_name = {'Baseline','B4','B6','B8','B11','B14','B17','B20','B23'};
% subj_name = {'subj6'};
% d = load_data(subj_name,block_name,folder,time);
% rotate = 0;
% data = analyze_data(d,block_name,0,rotate);

% time = 60; %in seconds
% block_name = {'B1_baseline','B3_dualBaseline','B4_darkBaseline','B5','B6','B7','B8','B9','B10','B11_dual','B12_dark','B13','B14','B15','B16','B17_dual','B18_dark','B19','B20','B21','B22','B23_dual','B24_dark','B25','B26','B27_dual','B28_dark'};
% graph_name = {'B1','B3 (D)','B4 (F)','B5','B6','B7','B8','B9','B10','B11 (D)','B12 (F)','B13','B14','B15','B16','B17 (D)','B18 (F)','B19','B20','B21','B22','B23 (D)','B24 (F)','B25','B26','B27 (D)','B28 (F)'};
% subj_name = {'subj7','subj8','subj9'};
% d = load_data(subj_name,block_name,folder,time);
% rotate = 0;
% data = analyze_data(d,block_name,0,rotate);

time = 60; %in seconds
block_name = {'B1_baseline','B3_dualBaseline','B4_darkBaseline','B5','B6','B7_dual','B8_dark','B9','B10','B11','B12','B13_dual','B14_dark','B15','B16','B17','B18','B19_dual','B20_dark','B22','B23','B24','B25_dual','B26_dark','B27','B28','B29','B30_dual','B31_dark','B32_habit','B33_habit'};
graph_name = {'B1','B3(D)','B4(F)','B5','B6','B7(D)','B8(F)','B9','B10','B11','B12','B13(D)','B14(F)','B15','B16','B17','B18','B19(D)','B20(F)','B22','B23','B24','B25(D)','B26(F)','B27','B28','B29','B30(D)','B31(F)','B32(H)','B33(H)'};
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
gblocks = [20];

% graph_amp_avg(data,block_name,gblocks,graph_name,'cursor'); % amplitude spectrums
% graph_MSE(data, block_name, graph_name); % mean squared error
% graph_lag(data, gblocks, graph_name,'Lhand'); % response lag
% graph_complexError(data, graph_name, gblocks); % complex tracking error
% graph_coherence(data, block_name, graph_name, 'cursor') % coherence
graph_xcorr(data,block_name,graph_name)