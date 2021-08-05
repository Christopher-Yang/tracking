function output = load_data(subj_name, block_name, folder, time, remove)
    
    disp('Loading...');
    most_freq = 0;
    samplingRate = 130;

    for i = 1:length(subj_name)
        disp(['   Subject ' num2str(i)]);
        for j = 1:length(block_name)
            if strcmp(subj_name{i},'subj14') && strcmp(block_name{j},'B24')
                path = [folder,subj_name{i},'/B23'];
            else
                path = [folder,subj_name{i},'/',block_name{j}];
            end
            tFile = dlmread([path,'/tFile.tgt']);
            Tb = 1/(tFile(1)/2);
            Nsamples = round(time*samplingRate);
            fnames = dir(path);
            full = [];
            
            for k = 1:length(fnames(not([fnames.isdir])))-1
                name = ['traj',num2str(k)];
                data.(name) = dlmread([path,'/',fnames(3+(k-1)).name],' ',6,0);
                idx = data.(name)(:,9)==1;
                data.(name) = data.(name)(idx,:);
                data.(name) = data.(name)(5*samplingRate+1:5*samplingRate+Nsamples,:);    %only trajectory after 5 sec warm up time is used
%                 data.(name) = cat(3,data.(name)(1:N,:),data.(name)(N+1:2*N,:),data.(name)(2*N+1:3*N,:)); % divide signal into trials with the length of base period
                full = cat(3,full,data.(name));
            end
            
            x = length(tFile)/6;
            if x > most_freq
                most_freq = x;
            end
            
%             bName = regexprep(block_name{j},'^...','');
            bName = block_name{j};
            
            if remove == 1
                if i == 4 && j == 1 % subj4, B1_baseline
                    output{i}.(bName).traj = full(:,:,2:5);
                else
                    output{i}.(bName).traj = full;
                end
            elseif remove == 2
                if i == 1 && j == 2 % 2-day, subj1, B2_darkBaseline
                    output{i}.(bName).traj = full(:,:,2:5);
                else
                    output{i}.(bName).traj = full;
                end
            else
                output{i}.(bName).traj = full;
            end
            output{i}.(bName).tFile = tFile;
        end
    end
%     output.most_freq = most_freq;
end