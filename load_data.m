% Extract raw data from the data files specified by the subject numbers
% (subj_rot and subj_mir), blocks (block_name), and file path (folder).
% This data is stored in "traj" while the frequencies, amplitudes, and
% phases used to generate target motion are stored in "tFile." See
% README.md for information what is contained in each data file and tFile.

function output = load_data(subj_rot, subj_mir, block_name, folder)

disp('Loading...');
time = 40; % seconds of data to be analyzed per trial; can be 20 or 40 secs
for p = 1:2 % iterate over groups of subjects
    if p == 1 % analyze rotation group
        subj_name = subj_rot;
    else % analyze mirror-reversal group
        subj_name = subj_mir;
    end
    for i = 1:length(subj_name) % iterate over subjects
        disp(['   ' subj_name{i}]);
        for j = 1:length(block_name) % iterate over blocks
            path = [folder,subj_name{i},'/',block_name{j}]; % set file path
            Nsamples = round(time*130.004); % number of samples to analyze
            fnames = dir(path); % find all filenames in path
            full = [];
            
            % read in frequencies, amplitudes, and phases of target
            tFile = dlmread([path,'/tFile.tgt']);
            
            % iterate over each file and put its contents into "full"
            for k = 1:length(fnames(not([fnames.isdir])))-1
                name = ['traj',num2str(k)];
                data.(name) = dlmread([path,'/',fnames(3+(k-1)).name] ...
                    ,' ',6,0);
                
                % only use final Nsamples worth of data to ensure clean
                % discrete Fourier transforms
                data.(name) = data.(name)(end-Nsamples+1:end,:);
                full = cat(3,full,data.(name));
            end
            
            % store loaded data
            if p == 1
                output.rot{i}.(block_name{j}).traj = full;
                output.rot{i}.(block_name{j}).tFile = tFile;
            else
                output.mir{i}.(block_name{j}).traj = full;
                output.mir{i}.(block_name{j}).tFile = tFile;
            end
        end
    end
end
end