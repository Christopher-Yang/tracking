function output = load_data(subj_name, folder, time)
    
    % display status 
    disp('Loading...');
    
    % fields from data file to be analyzed
    fields = {'tX_freq','tY_freq','tX_amp','tY_amp','tX_phase','tY_phase','cX_freq','cY_freq','cX_amp','cY_amp','cX_phase','cY_phase'};
    fields2 = {'time','cursorX','cursorY','targetX','targetY'};
    fields3 = {'extraTime','extraCursorX','extraCursorY','extraTargetX','extraTargetY'};
    
    allFields = [fields fields2 fields3];
    
    for i = 1:length(subj_name)
        disp(['   Subject ' num2str(i)]); % display progress
        path = [folder,subj_name{i}]; % file path for data
%         tFile = dlmread([path,'/tFile.csv']); 
%         fnames = dir(path);
        
        % read data
        d =  readtable([path '/testing6.csv']);
        
        % extract useful variables from the data table
        Ntrials = size(d,1);
        frameRate = d.frameRate(1);
        
        % set variables used in the analysis
        sampleIdx = 5*frameRate:(time+5)*frameRate-1; % data samples to be analyzed
%         Nsamples = length(sampleIdx); % number of samples to be analyzed
        tolerance = (1/frameRate) + 0.5/frameRate; % tolerance for defining duration of one frame
        count = 1; % index for storing sine parameters in cell array
        
        % organize data into appropriate format
        for k = 1:Ntrials
            trial = [];
            extraData = [];
            
            % format data from strings into matrices
            for j = 1:length(allFields)
                data = d.(allFields{j}){k}(2:end-1); % take of apostrophes and brackets from string
                data = strsplit(data,','); % separate data by commas
                data = cellfun(@str2double,data); % convert data into a matrix
                
                if sum(strcmp(allFields{j},fields))
                    sineParams.(allFields{j}){count} = data;
                elseif sum(strcmp(allFields{j},fields2))
                    trial = [trial data(sampleIdx)'];
                elseif sum(strcmp(allFields{j},fields3))
                    extraData = [extraData data(~isnan(data))'];
                end
            end
            
            % throw out extra data that was collected during 5-second ramp
            % period
            useSamples = find(extraData(:,1)>5);
            extraData = extraData(useSamples,:);

            % insert missing data (extraData) into the the main data
            % (trial)
            for m = 1:size(extraData,1)
                
                % find the time point immediately after and before missing
                % data
                after = find(trial(:,1)>=extraData(m,1),1);
                before = find(trial(:,1)<extraData(m,1),1,'last');
                
                replace = 1; % used to end while loop
                check = 1; % defines search breadth for NaNs
                while replace
                    % if NaN is first data point, place missing data there
                    if after == 1
                        trial(1,:) = extraData(m,:);
                        replace = 0;
                    
                    % if NaN is last data point, place missing data there
                    elseif before == size(trial,1)
                        trial(end,:) = extraData(m,:);
                        replace = 0;
                    
                    % if NaN is between before and after, place the missing
                    % data there
                    elseif after-before == 2
                        trial(before+1,:) = extraData(m,:);
                        replace = 0;
                    
                    % if there's 2 NaNs between before and after, fill
                    % missing data in index that's closer in time to before
                    % or after
                    elseif after-before == 3
                        if trial(after,1)-extraData(m,1) > extraData(m,1)-trial(before,1)
                            trial(before+1,:) = extraData(m,:);
                        else
                            trial(after-1,:) = extraData(m,:);
                        end
                        replace = 0;

                    % if more than 3 NaNs between before and after, don't
                    % fill data
                    elseif after-before > 3
                        disp('   Multiple NaNs between before and after')
                        replace = 0;
                        
                    % if NaN is outside of before and after, place the
                    % missing data as specified below
                    else
                        % find the indices of NaN(s)
                        idxRange = before-check:after+check; 
                        missing = idxRange(isnan(trial(idxRange,1)));
                        
                        % do different things based on number of NaNs found
                        switch length(missing)
                            % NaN not found
                            case 0
                                check = check + 1; % increase search range
                            
                            % one NaN found
                            case 1
                                % determine whether Nan is earlier or later
                                % in time
                                if missing < before
                                    idx = missing:before;
                                elseif missing > after
                                    idx = after:missing;
                                end
                                
                                % create "chunk," which will replace a
                                % subsection of the matrix
                                chunk = trial(idx,:);
                                chunk = sortrows([chunk; extraData(m,:)]);
                                chunk(end,:) = [];
                                trial(idx,:) = chunk;
                                
                                % end the while loop
                                replace = 0; 
                                
                            % two NaNs found
                            case 2
                                chunk = sortrows([trial(missing(1)+1:missing(2)-1,:); extraData(m,:)]);

                                % check whether to fill data by pushing
                                % data forward or backward
                                if ~isnan(trial(missing(1)-1,1)) % push data backward
                                    if chunk(1,1) - trial(missing(1)-1,1) < tolerance
                                        trial(missing(1):missing(2)-1,:) = chunk;
                                    end
                                elseif ~isnan(trial(missing(2)+1,1)) % push data forward
                                    if trial(missing(2)+1,1) - chunk(end,1) < tolerance
                                        trial(missing(1)+1:missing(2),:) = chunk;
                                    end
                                else % don't insert data if reasonable solution can't be found
                                    disp('Could not insert extra data');
                                end
                                
                                % end while loop
                                replace = 0;
                        end
                    end
                    if check > 100
                        disp(check);
                    end
                end
            end
            
            % check to see that data was inserted in chronological order
            sortCheck = trial(:,1);
            sortCheck = sortCheck(~isnan(sortCheck));
            if sum(sort(sortCheck) ~= sortCheck)
                error('Data insertion was done incorrectly');
            end
            
            % store data in output
            for j = 1:length(fields2)
                output{i}.(fields2{j})(:,k) = trial(:,j);
            end
            
            count = count + 1;
        end
        
        frameDrops = sum(isnan(output{i}.time),1)./frameRate; % calculate frame drops for each trial
        
        firstTime = nanmean(output{i}.time(1,:));
        for j = 1:Ntrials
            if ~isnan(output{i}.time(1,j))
                output{i}.time(:,j) = output{i}.time(:,j) - output{i}.time(1,j); % begin time at 0
            else
                output{i}.time(:,j) = output{i}.time(:,j) - firstTime;
            end
        end
        
        for j = 1:length(fields)
            output{i}.(fields{j}) = sineParams.(fields{j});
        end
        output{i}.trialType = d.trialType;
        output{i}.frameRate = frameRate;
        output{i}.frameDrops = frameDrops;
        output{i}.OS = d.OS{1};
        output{i}.browser = d.browser{1};
        output{i}.cmConvert = d.cmConvert(1);
    end
end

