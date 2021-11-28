% This file is associated with the article "Emergence of habitual control
% in a novel motor skill over multiple days of practice". Running
% main.m generates plots for Figures 3, 5, and S1.
% 
% Part 1 loads raw data from the tracking task and performs initial data 
% analysis, storing the data in a structure called "data". This
% structure is organized as follows: data.(group){subject}.(block). See
% analyze_data() for more details on subfields in "data". Part 2 performs 
% follow-up data analysis and plots figures.
% 
% All analyses can be performed by running main.m

%% PART 1: INITIAL DATA  ANALYSIS
clear all;

% analysis for 2-day group
folder = 'Data/denovo_2day/';
block_name.day2 = {'B1_baseline','B2_darkBaseline','B3','B6_dark','B9',...
    'B10_dark','B11_habit','B12_habit','B13_habitDark'};
blockType.day2 = [1 2 1 2 1 2 3 3 2];

% graph_name.day2 = {'Baseline','Baseline (D)','Early','Day 1 (D)',...
%     'Day 2','Day 2 (D)','Flip1 (F)','Flip2 (F)','Flip3 (FD)'};
subj_name = {'subj1','subj3','subj4','subj5','subj6','subj7','subj8',...
    'subj9','subj10','subj11','subj12','subj13','subj14'};
remove = 2; % flag to remove bad data from subj1 B2_darkBaseline
d = load_data(subj_name,block_name.day2,folder,remove);
data.day2 = analyze_data(d,block_name.day2);

% analysis for 5-day group
folder = 'Data/denovo_5day/';
block_name.day5 = block_name.day2(1:end-5);
blockType.day5 = blockType.day2(1:end-3);
block_name.day5 = [block_name.day5 {'B10','B11_dark','B15','B16_dark',...
    'B20','B21_dark','B24','B25_dark','B26_habit','B27_habit',...
    'B28_habitDark'}];
blockType.day5 = [blockType.day5 1 2 1 2 1 2 3 3 2];
% graph_name.day5 = [graph_name.day5 {'Day 3','Day 3 (D)','Day 4',...
%     'Day 4 (D)','Day 5','Day 5 (D)','Flip (F)','Flip2 (F)','Flip3 (FD)'}];
subj_name = {'subj13','subj15','subj17','subj18','subj19','subj21',...
    'subj22','subj23','subj24','subj25','subj26','subj27','subj28',...
    'subj29'};
remove = 0; % don't need to remove any data from this group
d = load_data(subj_name,block_name.day5,folder,remove);
data.day5 = analyze_data(d,block_name.day5);

% analysis for 10-day group
folder = 'Data/denovo_10day/';
block_name.day10 = block_name.day5(1:end-5);
blockType.day10 = blockType.day5(1:end-3);
block_name.day10 = [block_name.day10 {'B25','B26_dark','B30','B31_dark',...
    'B35','B36_dark','B40','B41_dark','B45','B46_dark','B49','B50_dark',...
    'B51_habit','B52_habit','B53_habitDark'}];
blockType.day10 = [blockType.day10 1 2 1 2 1 2 1 2 1 2 3 3 2];
% graph_name.day10 = [graph_name.day10 {'Day 6','Day 6 (D)','Day 7',...
%     'Day 7 (D)','Day 8','Day 8 (D)','Day 9','Day 9 (D)','Day 10',...
%     'Day 10 (D)','Flip (F)','Flip2 (F)','Flip3 (FD)'}];
subj_name = {'subj1','subj2','subj3','subj4','subj5'};
remove = 1; % flag to remove bad data from subj4 B1_baseline
d = load_data(subj_name,block_name.day10,folder,remove);
data.day10 = analyze_data(d,block_name.day10);

disp('Done')

%% PART 2: FOLLOW-UP DATA ANALYSIS AND PLOTTING FIGURES

% Figure 3A
graph_cursor(data)

% Figure 3B
graph_coherence(data, block_name, blockType)

% Figures 3C and 5A
graph_gainMatrix(data, block_name, blockType);

% Figure 5B
graph_flipGain(data)