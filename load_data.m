function output = load_data(subj_rot, subj_rot_i, block_name, folder)
    
    disp('Loading...');
    time = 40; % seconds of data to be analyzed per trial; can be 20 or 40 secs
    for p = 1:2
        if p == 1 % analyze rotation group
            subj_name = subj_rot;
        else % analyze mirror-reversal group
            subj_name = subj_rot_i;
        end
        for i = 1:length(subj_name)
            disp(['   Subject ' num2str(i)]);
            for j = 1:length(block_name)
                path = [folder,subj_name{i},'/',block_name{j}]; % set file path
                tFile = dlmread([path,'/tFile.tgt']); % read in frequencies, amplitudes, and phases of target
                Nsamples = round(time*130.004); % number of samples to analyze
                fnames = dir(path); % find all filenames in path
                full = []; 
                
                for k = 1:length(fnames(not([fnames.isdir])))-1
                    name = ['traj',num2str(k)];
                    data.(name) = dlmread([path,'/',fnames(3+(k-1)).name],' ',6,0); % read in all data
                    data.(name) = data.(name)(end-Nsamples+1:end,:);    % only use final Nsamples worth of data
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