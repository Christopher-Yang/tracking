function [data,d,groups,block_name,graph_name,theta] = loadSubj(x)

switch x
    
    case 1 % VMR vs MR with 150 P2P
        folder = 'Data/VMR90_vs_MR/';
        time = 40; %in seconds
        theta = 90;
        groups = {'rot','rot_i'};
        block_name = {'no_rot1','rot1','rot2','rot3','rot4','no_rot2'};
        graph_name = {'Baseline','Early','Train1','Train2','Late','Post'};
        subj_name = {'subj14','subj15','subj16','subj17','subj18','subj19','subj21','subj22','subj23','subj24','subj25','subj26','subj27','subj28','subj29','subj30','subj31','subj32','subj33','subj34'};
        subj_rot = {'subj17','subj18','subj21','subj22','subj24','subj25','subj28','subj31','subj32','subj33'};
        subj_rot_i = {'subj14','subj15','subj16','subj19','subj23','subj26','subj27','subj29','subj30','subj34'};
        d = load_data(subj_name,block_name,folder,time);
        data.rot = analyze_data(d,subj_rot,block_name,false,1);
        data.rot_i = analyze_data(d,subj_rot_i,block_name,false,2);

    case 2 % VMR vs MR with 10 P2P 
        folder = 'Data/10P2P/';
        time = 40; %in seconds
        theta = 90;
        groups = {'rot','rot_i'};
        block_name = {'no_rot1','rot1','rot2','rot3','rot4','no_rot2'};
        graph_name = {'Baseline','Early','Train1','Train2','Late','Post'};
        subj_name = {'subj37','subj38','subj39','subj40','subj41','subj42','subj43','subj44'};
        subj_rot = {'subj37','subj38','subj39','subj40','subj41'};
        subj_rot_i = {'subj42','subj43','subj44'};
        d = load_data(subj_name,block_name,folder,time);
        data.rot = analyze_data(d,subj_rot,block_name,false,3);
        data.rot_i = analyze_data(d,subj_rot_i,block_name,false,4);
        
    case 3 % control2- testing why x- and y-hand responses are asymmetrical
        folder = 'Data/XY_asymmetry/';
        time = 40; %in seconds
        theta = 90;
        groups = {'rot','rot_i'};
        block_name1 = {'XY_laptop1','XY_laptop2','XY_desktop1','XY_desktop2'};
        block_name2 = {'XY_rot1','XY_rot2'};
        graph_name = {'laptop 1','laptop 2','desktop 1','desktop 2','diagonal 1','diagonal 2'};
        subj_name = {'subj45','subj46','subj47','subj48'};
        subj_rot = {'subj45','subj46','subj47','subj48'};
        subj_rot_i = {'subj45','subj46','subj47','subj48'};
        d1 = load_data(subj_name,block_name1,folder,time);
        d2 = load_data(subj_name,block_name2,folder,time);
        data.rot = analyze_data(d1,subj_rot,block_name1,false);
        data.rot_i = analyze_data(d2,subj_rot_i,block_name2,false);
        
    case 4 % control3 - my data for pure X, pure Y, and new rotated sum of sines
        folder = 'Data/pilot/';
        time = 40; %in seconds
        theta = 90;
        groups = {'rot'};
        block_name = {'rot'};
        graph_name = {'rot'};
        subj_name = {'test'};
        subj_rot = {'test'};
        d = load_data(subj_name,block_name,folder,time);
        data.rot = analyze_data(d,subj_rot,block_name,true);
        
    case 5 % VMR vs MR w/ no P2P training
        folder = 'Data/no_P2P/';
        time = 40; %in seconds
        theta = 90;
        groups = {'rot','rot_i'};
        block_name = {'baseline','rot1','rot2','rot3','rot4','after'};
        graph_name = {'Baseline','Early','Train1','Train2','Late','Aftereffects'};
        subj_name = {'subj1','subj2','subj3','subj4','subj5','subj6','subj7','subj8','subj9','subj10','subj11','subj12','subj13','subj14','subj15','subj16','subj17','subj18','subj19','subj20'};
        subj_rot = {'subj2','subj5','subj7','subj8','subj9','subj10','subj13','subj15','subj17','subj19'};
        subj_rot_i = {'subj1','subj3','subj4','subj6','subj11','subj12','subj14','subj16','subj18','subj20'};
        d = load_data(subj_name,block_name,folder,time);
        data.rot = analyze_data(d,subj_rot,block_name,false,0);
        data.rot_i = analyze_data(d,subj_rot_i,block_name,false,0);
        
    case 6 % VMR15 with no P2P training
        folder = 'Data/no_P2P_15deg/';
        time = 40; %in seconds
        theta = 15;
        groups = {'rot'};
        block_name = {'baseline','rot1','rot2','rot3','rot4','after'};
        graph_name = {'Baseline','Early','rot2','rot3','Late','Post'};
        % block_name = {'a'}; % I used different sinusoids here with less spectral density
        % block_name = {'test'};
%         subj_name = {'subj1','subj2','subj3','subj4','subj5','subj6','subj7','subj8','subj9','subj10'};
        subj_name = {'subj11','subj12','subj13','subj14','subj15','subj17','subj18','subj20','subj21','subj22'};
        theta = 15;
        subj_rot = subj_name;
        d = load_data(subj_name,block_name,folder,time);
        data.rot = analyze_data(d,subj_rot,block_name,false,0);
        
    case 7 % MR with 15 deg VMR added after 1 day of training
        folder = 'Data/mir_smallrot/';
        time = 40; %in seconds
        theta = 90;
        groups = {'rot_i'};
        block_name = {'D1baseline','D1rot1','D1rot2','D1rot3','D1rot4','D1adapt','D1after'};
        graph_name = {'D1baseline','D1rot1','D1rot2','D1rot3','D1rot4','D1adapt','D1after'};
        subj_name = {'subj4','subj6','subj7','subj8'};
        subj_rot = subj_name;
        theta = 90;
        d = load_data(subj_name,block_name,folder,time);
        data.rot_i = analyze_data(d,subj_rot,block_name,false,0);

    case 8 % VMR90 with no cursor feedback during tracking but with feedback during 150 P2P
        folder = 'Data/dark/';
        time = 40; %in seconds
        theta = 90;
        groups = {'rot'};
        block_name = {'baseline','dark_baseline','rot1','rot2','rot3','rot4','after'};
        graph_name = {'baseline','dark_baseline','rot1','rot2','rot3','rot4','after'};
        subj_name = {'subj1','subj2','subj3','subj4','subj5'};
        subj_rot = subj_name;
        d = load_data(subj_name,block_name,folder,time);
        data.rot = analyze_data(d,subj_rot,block_name,false,0);
        
    case 9 % dual task pilot
        groups = {'rot'};
        block_name = {'avg','avg_dual'};
        graph_name = {'avg','avg dual'};
        subj_name = {'subj1'};
        subj_rot = {'subj1'};
        d = load_data(subj_name,block_name,folder,time);
        data.rot = analyze_data(d,subj_rot,block_name,false,0);
        
    case 10
        folder = 'Data/no_P2P_15deg/';
        time = 40; %in seconds
        theta = 15;
        groups = {'rot'};
        block_name = {'baseline','rot1','rot2','rot3','rot4','after'};
        graph_name = {'Baseline','Early','rot2','rot3','Late','Post'};
        subj_name = {'subj11','subj12','subj13','subj14','subj15','subj17','subj18','subj19','subj20','subj21'};
        theta = 15;
        subj_rot = subj_name;
        d = load_data(subj_name,block_name,folder,time);
        data.rot = analyze_data(d,subj_rot,block_name,false,0);
end