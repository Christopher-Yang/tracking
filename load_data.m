function output = load_data(subj_rot, subj_rot_i, block_name, folder)
    
    disp('Loading...');
    time = 40; % seconds of data to be analyzed per trial; can be 20 or 40 secs
    for p = 1:2
        if p == 1
            subj_name = subj_rot;
        else
            subj_name = subj_rot_i;
        end
        for i = 1:length(subj_name)
            disp(['   Subject ' num2str(i)]);
            for j = 1:length(block_name)
                path = [folder,subj_name{i},'/',block_name{j}];
                tFile = dlmread([path,'/tFile.tgt']);
                Tb = 1/(tFile(1)/2);
                Ncycles = floor(time/Tb);
                trial_time = Tb*Ncycles;
                Nsamples = round(trial_time*130.004);
                fnames = dir(path);
                full = [];
                
                for k = 1:length(fnames(not([fnames.isdir])))-1
                    name = ['traj',num2str(k)];
                    data.(name) = dlmread([path,'/',fnames(3+(k-1)).name],' ',6,0);
                    data.(name) = data.(name)(end-Nsamples+1:end,:);    %only trajectory after 5 sec warm up time is used
                    full = cat(3,full,data.(name));
                end
                
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