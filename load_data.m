function output = load_data(subj_name, block_name, folder,time)
    
    disp('Loading...');
    most_freq = 0;

    for i = 1:length(subj_name)
        disp(['   ' subj_name{i}]);
        for j = 1:length(block_name)
            path = [folder,subj_name{i},'/',block_name{j}];
            tFile = dlmread([path,'/tFile.tgt']);
            Tb = 1/(tFile(1)/2);
%             Ncycles = floor(time/Tb);
%             trial_time = Tb*Ncycles;
%             Nsamples = round(trial_time*130.004)+13;
            Nsamples = round(time*130.004)+13;
            fnames = dir(path);
            full = [];
            
            for k = 1:length(fnames(not([fnames.isdir])))-1
                name = ['traj',num2str(k)];
                data.(name) = dlmread([path,'/',fnames(3+(k-1)).name],' ',6,0);   
                data.(name) = data.(name)(end-Nsamples:end,:);    %only trajectory after 5 sec warm up time is used
                full = cat(3,full,data.(name));
            end
            
            x = length(tFile)/6;
            if x > most_freq
                most_freq = x;
            end
            
            output.(subj_name{i}).(block_name{j}).traj = full;
            output.(subj_name{i}).(block_name{j}).tFile = tFile;
        end
    end
    output.most_freq = most_freq;
end