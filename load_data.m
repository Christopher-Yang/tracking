% Extract raw data from the data files specified by the subject numbers
% (subj_rot and subj_mir), blocks (block_name), and file path (folder).
% This data is stored in "traj" while the frequencies, amplitudes, and
% phases used to generate target motion are stored in "tFile." See
% README.md for information what is contained in each data file and tFile.

function output = load_data(subj_rot, subj_mir, block_name, folder)

disp('Loading...');
time = 40; % seconds of data to be analyzed per trial; can be 20 or 40 secs
groups = {'rot','mir'};
for p = 1
% for p = 1:2 % iterate over groups of subjects
    if p == 1 % analyze rotation group
        subj_name = subj_rot;
    else % analyze mirror-reversal group
        subj_name = subj_mir;
    end
    for i = 1:length(subj_name) % iterate over subjects
        disp(['   ' subj_name{i}]);
        for j = 1:length(block_name) % iterate over blocks
            
            path = [folder,subj_name{i},'/',block_name{j}]; % set file path
            Nsamples = round(time*130); % number of samples to analyze
            fnames = dir(path); % find all filenames in path
            full = [];
            
            % read in frequencies, amplitudes, and phases of target
            tFile = dlmread([path,'/tFile.tgt']);
            
            % iterate over each file and put its contents into "full"
            Ntrials = length(fnames(not([fnames.isdir])))-1;
            trialType = NaN(Ntrials,1);
            for k = 1:Ntrials
                name = ['traj',num2str(k)];
                data.(name) = dlmread([path,'/',fnames(3+(k-1)).name] ...
                    ,' ',10,0);
                
                % only use final Nsamples worth of data to ensure clean
                % discrete Fourier transforms
                data.(name) = data.(name)(end-Nsamples+1:end,:);
                full = cat(3,full,data.(name));
                
                % open the file
                fid=fopen([path '/' fnames(3+(k-1)).name]);
                
                % set linenum to the desired line number that you want to import
                linenum = 4;
                
                % use '%s' if you want to read in the entire line or use '%f' if you want to read only the first numeric value
                line4 = textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum-1);
                line5 = textscan(fid,'%s',1,'delimiter','\n');
                line7 = textscan(fid,'%s',1,'delimiter','\n', 'headerlines',1);
                fclose(fid);
                
                idx = isstrprop(line4{1}{1},'digit');
                trialType(k) = str2double(line4{1}{1}(idx));
            end
            
            idx = isstrprop(line5{1}{1},'digit');
            bimanual_mode = str2double(line5{1}{1}(idx));
            
            idx = isstrprop(line7{1}{1},'digit');
            rotation = str2double(line7{1}{1}(idx));
            
            % store loaded data
            output.(groups{p}){i}.(block_name{j}).traj = full;
            output.(groups{p}){i}.(block_name{j}).tFile = tFile;
            output.(groups{p}){i}.(block_name{j}).trialType = trialType;
            output.(groups{p}){i}.(block_name{j}).bimanual_mode = bimanual_mode;
            output.(groups{p}){i}.(block_name{j}).rotation = rotation;
        end
    end
end
end