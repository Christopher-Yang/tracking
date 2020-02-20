clear all;

folder = 'Data/denovo_5day/';

% time = 40; %in seconds
% block_name = {'B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11','B12','B13','B14','B15','B16','B17','B18','B19','B21','B22_habit','B23','B24'};
% graph_name = {'B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11','B12','B13','B14','B15','B16','B17','B18','B19','B21','B22 (H)','B23','B24'};
% block_name = {'B1','B2','B8','B14','B17','B18','B19','B21','B22_habit'};
% graph_name = {'B1','B2','B8','B14','B17','B18','B19','B21','B22 (H)'};
% subj_name = {'subj1','subj2','subj3'};
% d = load_data(subj_name,block_name,folder,time);
% rotate = 0;
% data = analyze_data(d,block_name,0,rotate);

% subj 4 had trial lengths too short; subj 5 was poor at task
% time = 60; %in seconds
% block_name = {'B1_baseline','B4','B11','B17','B23','B27','B31_habit','B32_habit'};
% graph_name = {'B1','B4','B11','B17','B23','B27','B31 (H)','B32 (H)'};
% block_name = {'B1_baseline','B3_dualBaseline','B6','B7_dual','B11','B12_dual','B17','B18_dual','B23','B24_dual','B27','B28_dual'};
% graph_name = {'B1','B3 (D)','B6','B7 (D)','B11','B12(D)','B17','B18 (D)','B23','B24 (D)','B27','B28 (D)'};
% block_name = {'B1_baseline','B4','B11','B13_dark','B17','B19_dark','B23','B25_dark','B27','B29_dark'};
% graph_name = {'B1','B4','B11','B13 (F)','B17','B19 (F)','B23','B25 (F)','B27','B29 (F)'};
% block_name = {'B1_baseline','B3_dualBaseline','B11','B12_dual','B17','B18_dual','B23','B24_dual','B27','B28_dual'};
% graph_name = {'B1','B3 (D)','B11','B12 (D)','B17','B18 (D)','B23','B24 (D)','B27','B28 (D)'};
% subj_name = {'subj6'};
% d = load_data(subj_name,block_name,folder,time);
% rotate = 0;
% data = analyze_data(d,block_name,0,rotate);

% for these subjects, I rotated the bimanual mapping so outward and up made
% the cursor go inward and up but the target frequencies weren't aligned to
% the hand motion
% time = 60; %in seconds
% block_name = {'B1_baseline','B3_dualBaseline','B4_darkBaseline','B5','B6','B7','B8','B9','B10','B11_dual','B12_dark','B13','B14','B15','B16','B17_dual','B18_dark','B19','B20','B21','B22','B23_dual','B24_dark','B25','B26','B27_dual','B28_dark'};
% graph_name = {'B1','B3 (D)','B4 (F)','B5','B6','B7','B8','B9','B10','B11 (D)','B12 (F)','B13','B14','B15','B16','B17 (D)','B18 (F)','B19','B20','B21','B22','B23 (D)','B24 (F)','B25','B26','B27 (D)','B28 (F)'};
% block_name = {'B1_baseline','B5','B10','B22','B26'};
% graph_name = {'B1','B5','B19','B22','B26'};
% subj_name = {'subj7','subj8','subj9'};
% block_name = {'B1_baseline','B5','B10','B22','B26','B29_habit','B30_habit'};
% graph_name = {'B1','B5','B19','B22','B26','B29 (H)','B30 (H)'};
% subj_name = {'subj8','subj9'};
% d = load_data(subj_name,block_name,folder,time);
% rotate = 0;
% data = analyze_data(d,block_name,0,rotate);

% these received the proper target motion for the bimanual mapping used for
% subj7-9
% time = 60; %in seconds
% block_name = {'B1_baseline','B3_dualBaseline','B4_darkBaseline','B5','B6','B7_dual','B8_dark','B9','B10','B11','B12','B13_dual','B14_dark','B15','B16','B17','B18','B19_dual','B20_dark','B22','B23','B24','B25_dual','B26_dark','B27','B28','B29','B30_dual','B31_dark','B32_habit','B33_habit'};
% graph_name = {'B1','B3(D)','B4(F)','B5','B6','B7(D)','B8(F)','B9','B10','B11','B12','B13(D)','B14(F)','B15','B16','B17','B18','B19(D)','B20(F)','B22','B23','B24','B25(D)','B26(F)','B27','B28','B29','B30(D)','B31(F)','B32(H)','B33(H)'};
% block_name = {'B1_baseline','B4_darkBaseline','B6','B8_dark','B12','B14_dark','B18','B20_dark','B24','B26_dark','B29','B31_dark'};
% graph_name = {'B1','B4 (F)','B6','B8 (F)','B12','B14 (F)','B18','B20 (F)','B24','B26 (F)','B29','B31 (F)'};
% block_name = {'B1_baseline','B5','B12','B24','B29','B32_habit','B33_habit'};
% graph_name = {'B1','B5','B12','B24','B29','B32 (H)','B33 (H)'};
% block_name = {'B1_baseline','B3_dualBaseline','B6','B7_dual','B12','B13_dual','B18','B19_dual','B24','B25_dual','B29','B30_dual'};
% graph_name = {'B1','B3 (D)','B6','B7 (D)','B12','B13 (D)','B18','B19 (D)','B24','B25 (D)','B29','B30 (D)'};
% subj_name = {'subj10','subj11'};
% d = load_data(subj_name,block_name,folder,time);
% rotate = 1;
% data = analyze_data(d,block_name,0,rotate);

% received the normal bimanual mapping
% time = 60; %in seconds
% % block_name = {'B1_baseline','B2_dark','B3','B4','B5','B6_dark','B7','B8','B9','B10','B11_dark','B12','B13','B14','B15','B16_dark','B17','B18','B19','B20','B21_dark','B22','B23','B24','B25_dark','B26_habit','B27_habit'};
% % graph_name = {'B1','B2 (D)','B3','B4','B5','B6 (D)','B7','B8','B9','B10','B11 (D)','B12','B13','B14','B15','B16(D)','B17','B18','B19','B20','B21 (D)','B22','B23','B24','B25 (D)','B26 (H)','B27 (H)'};
% block_name = {'B1_baseline','B3','B10','B15','B20','B24','B26_habit','B27_habit'};
% graph_name = {'Baseline','Early','Day 2','Day 3','Day 4','Day 5','Flip 1 (F)','Flip 2 (F)'};
% % block_name = {'B1_baseline','B3','B10','B20','B23','B26_habit','B27_habit','B28_habitDark'};
% % graph_name = {'B1','B3','B10','B20','B23','B26 (H)','B27 (H)','B27 (HD)'};
% % block_name = {'B1_baseline','B2_dark','B10','B11_dark','B15','B16_dark','B20','B21_dark','B24','B25_dark'};
% % graph_name = {'B1','B2 (F)','B10','B11 (F)','B15','B16 (F)','B20','B21 (F)','B24','B25 (F)'};
% subj_name = {'subj13','subj14','subj15','subj16','subj17','subj18','subj19','subj20','subj21','subj22'};
% % subj_name = {'subj13','subj15','subj16','subj17','subj18','subj19','subj20'};
% d = load_data(subj_name,block_name,folder,time);
% rotate = 0;
% data = analyze_data(d,block_name,0,rotate);

folder = 'Data/denovo_10day/';
time = 60; %in seconds
% block_name = {'B1_baseline','B3','B20','B35','B49','B51_habit','B52_habit','B53_habitDark'};
% graph_name = {'Baseline','Early','Day 4','Day 7','Day 10','Flip 1','Flip 2','Flip (dark)'};
% block_name = {'B6_dark','B11_dark','B16_dark','B21_dark','B26_dark','B31_dark','B36_dark','B41_dark','B46_dark','B50_dark','B53_habitDark'};
% graph_name = {'B6','B11','B16','B21','B26','B31','B36','B41','B46','B50','Flip (dark)'};
block_name = {'B1_baseline','B3','B10','B15','B20','B25','B30','B35','B40','B45','B49'};
graph_name = {'Baseline','Early','Day 2','Day 3','Day 4','Day 5','Day 6','Day 7','Day 8','Day 9','Day 10'};
subj_name = {'subj1','subj2','subj3','subj4'};
remove = 1; % subj4 B1_baseline the Flock of Birds was bad
d = load_data(subj_name,block_name,folder,time,remove);
rotate = 0;
data = analyze_data(d,block_name,0,rotate);

% folder = 'Data/denovo_2day/';
% time = 60; %in seconds
% block_name = {'B1_baseline','B3','B5','B9','B11_habit','B12_habit'};
% graph_name = {'Baseline','Early','Day 1','Day 2','Flip 1 (F)','Flip 2 (F)'};
% % block_name = {'B2_darkBaseline','B6_dark','B10_dark'};
% % graph_name = {'Baseline','Day 1','Day 2'};
% subj_name = {'subj1','subj2','subj3','subj4','subj5','subj6','subj7'};
% d = load_data(subj_name,block_name,folder,time);
% rotate = 0;
% data = analyze_data(d,block_name,0,rotate);

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