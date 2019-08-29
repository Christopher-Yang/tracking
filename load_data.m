function output = load_data(subj_name, block_name, folder,time)
    
    disp('Loading...');
    most_freq = 0;

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
            N = Nsamples/Ncycles;
            
            for k = 1:length(fnames(not([fnames.isdir])))-1
                name = ['traj',num2str(k)];
                data.(name) = dlmread([path,'/',fnames(3+(k-1)).name],' ',6,0);
                idx = data.(name)(:,9)==1;
                data.(name) = data.(name)(idx,:);
                data.(name) = data.(name)(end-Nsamples+1:end,:);    %only trajectory after 5 sec warm up time is used
%                 data.(name) = cat(3,data.(name)(1:N,:),data.(name)(N+1:2*N,:),data.(name)(2*N+1:3*N,:)); % divide signal into trials with the length of base period
                full = cat(3,full,data.(name));
            end
            
            x = length(tFile)/6;
            if x > most_freq
                most_freq = x;
            end
            
%             bName = regexprep(block_name{j},'^...','');
            bName = block_name{j};
            
            output{i}.(bName).traj = full;
            output{i}.(bName).tFile = tFile;
        end
    end
%     output.most_freq = most_freq;
end