function output = load_data(subj_name, folder, time, remove)
    
    disp('Loading...');
    most_freq = 0;
    samplingRate = 60;
    fields = {'time','cursorX','cursorY','targetX','targetY','cX_input','cY_input','tX_input','tY_input'};
    fields2 = {'time','cursorX','cursorY','targetX','targetY'};

    for i = 1:length(subj_name)
        disp(['   Subject ' num2str(i)]);
        path = [folder,subj_name{i}];
        tFile = dlmread([path,'/tFile.csv']);
        Nsamples = round(time*samplingRate);
        fnames = dir(path);
        full = [];
        
        d =  readtable([path '/firefox_multipletrials2.csv']);
        Ntrials = size(d,1);
        
        for j = 1:length(fields)
            field = d.(fields{j});
            for k = 1:Ntrials
                trial = field{k}(2:end-1);
                trial = strsplit(trial,',');
                trial = cellfun(@str2double,trial);
                if sum(strcmp(fields2,fields{j}))
                    output{i}.(fields{j})(:,k) = trial(end-Nsamples+1:end);
                else
                    output{i}.(fields{j})(:,k) = trial;
                end
            end
        end
        
        output{i}.cursorPerturb = d.cursorSines;
        
%         for j = 1:length(fields2)
%             output{i}.(fields2{j}) = output{i}.(fields2{j})(end-Nsamples+1:end);
%         end
        
        output{i}.time = output{i}.time - repmat(output{i}.time(1,:),[Nsamples 1]); % begin time at 0
        timeDiff = diff(output{i}.time,1); % calculate period between samples
        idx = timeDiff > 1/samplingRate + 0.01; % find time points where frames were dropped
        frameDrops = timeDiff.*idx - (1/samplingRate).*idx; % subtract off normal sampling period from frame drop
        frameDrops = sum(frameDrops,1); % sum frame drops for each trial
        
        output{i}.frameDrops = frameDrops;
        
    end
end

% function mat = convertMatrix(cells)
% Ntrials = size(cells,1);
% for k = 1:Ntrials
%     trial = cells{k}(2:end-1);
%     trial = strsplit(trial,',');
%     mat(:,k) = cellfun(@str2double,trial);
% %     if strcmp(fields{j},'time') || strcmp(fields{j},'cursorX') || strcmp(fields{j},'cursorY') || strcmp(fields{j},'targetX') || strcmp(fields{j},'targetY')
% %         output{i}.(fields{j})(:,k) = trial(end-Nsamples+1:end);
% %     end
% end
% end