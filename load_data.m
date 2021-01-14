function output = load_data(folder, time)
    
    % display status 
    disp('Loading...');
    
    % fields from data file to be analyzed
    fields = {'tX_freq','tY_freq','tX_amp','tY_amp','tX_phase','tY_phase','cX_freq','cY_freq','cX_amp','cY_amp','cX_phase','cY_phase'};
    fields2 = {'time','cursorX','cursorY','handX','handY','targetX','targetY','cursorX_input','cursorY_input'};
    fields3 = {'extraTime','extraCursorX','extraCursorY','extraHandX','extraHandY','extraTargetX','extraTargetY','extraCursorX_input','extraCursorY_input'};

    allFields = [fields fields2 fields3];
    fnames = dir(folder);
    Nsubj = length(fnames)-2;
    
    for i = 1:Nsubj
        disp(['   Subject ' num2str(i)]); % display progress
        
        % read data
        d = readtable([folder '/' fnames(2+i).name]);
        
        % only analyze tracking trials
        tracking = strcmp(d.task,'tracking');
        d = d(tracking,:);
        
        % extract useful variables from the data table
        Ntrials = size(d,1);
        frameRate = d.frameRate(1);
        
        % set variables used in the analysis
        sampleIdx = 5*frameRate:(time+5)*frameRate-1; % data samples to be analyzed
        Nsamples = length(sampleIdx); % number of samples to be analyzed
        tolerance = (1/frameRate) + 0.5/frameRate; % tolerance for defining duration of one frame
        count = 1; % index for storing sine parameters in cell array
        count2 = 0; % index for arranging data into blocks
        block = 1;
        
        % organize data into appropriate format
        for k = 1:Ntrials
            trial = [];
            extraData = [];
            
            % increments count2 if block increases
            if block == d.block(k) % current block is still the same
                count2 = count2 + 1;
            else % current block has incremented
                block = d.block(k);
                count2 = 1;
            end
            
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
            if ~isempty(extraData)
                useSamples = find(extraData(:,1)>5);
                extraData = extraData(useSamples,:);
            end
            
            % check whether there's any missing data
            Nmissing = sum(sum(isnan(trial)))/size(trial,2);
            Nextra = size(extraData,1);
            
            % insert missing data (extraData) into the the main data
            % (trial)
            if Nmissing ~= 0 % skip insertion if there's no missing data
                if Nextra <= Nmissing % this handles insertion when number of extra data <= missing data
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
                                disp('      Multiple NaNs between before and after')
                                replace = 0;
                                
                            % if NaN is outside of before and after, place the
                            % missing data as specified below
                            else
                                % find the indices of NaN(s)
                                if before-check < 1
                                    first = 1;
                                else
                                    first = before-check;
                                end
                                if after+check > Nsamples
                                    last = Nsamples;
                                else
                                    last = after+check;
                                end
                                idxRange = first:last;
                                
                                try
                                    missing = idxRange(isnan(trial(idxRange,1)));
                                catch
                                    error('      Something went wrong')
                                end
                                
                                % do different things based on number of NaNs found
                                switch length(missing)
                                    
                                    case 0 % NaN not found
                                        if check > 18 % only search within 300 ms of desired data position
                                            disp('      NaN could not be found');
                                            replace = 0;
                                        else % increase search range
                                            check = check + 1;
                                        end
                                        
                                        
                                    case 1 % one NaN found
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
                                        
                                        
                                    case 2 % two NaNs found
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
                                            disp('      Could not find reasonable place to insert data');
                                        end
                                        
                                        % end while loop
                                        replace = 0;
                                end
                            end
                        end
                    end
                else % this best handles filling in missing data if frame rate >> 60 Hz
                    % find NaNs in data
                    idx2 = find(isnan(trial(:,1)));
                    
                    % loop over all NaNs
                    for m = 1:length(idx2)
                        
                        % search backward for first non-NaN value
                        before = idx2(m)-1;
                        search = true;
                        while search
                            if isnan(trial(before))
                                before = before - 1;
                            else
                                search = false;
                            end
                        end
                        
                        % search forward for first non-NaN value
                        after = idx2(m)+1;
                        search = true;
                        while search
                            if isnan(trial(after))
                                after = after + 1;
                            else
                                search = false;
                            end
                        end
                        
                        Nidx = after - before + 1; % number of indices between after and before + 2
                        fillIdx = idx2(m) - before + 1; % the index of interest
                        
                        % linearly interpolate between before and after
                        sample_interp = linspace(trial(before), trial(after), Nidx);
                        
                        % find which extra data most closely matches
                        % interpolated data
                        [~, idx3] = min(abs(extraData(:,1) - sample_interp(fillIdx)));
                        
                        % if extra data that's chosen falls between the
                        % times for "before" and "after," then fill in
                        % missing data
                        if extraData(idx3,1) < trial(after,1) && extraData(idx3,1) > trial(before,1)
                            trial(idx2(m),:) = extraData(idx3,:);
                        end
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
                output{i}{d.block(k)}.(fields2{j})(:,count2) = trial(:,j);
            end
            
            count = count + 1;
        end
        
        Nblock = max(d.block);
        
        for k = 1:max(d.block)
            frameDrops = sum(isnan(output{i}{k}.time),1)./frameRate; % calculate frame drops for each trial
            
            % adjust time so times start at 0
            firstTime = nanmean(output{i}{k}.time(1,:));
            for j = 1:Ntrials/Nblock
                if ~isnan(output{i}{k}.time(1,j))
                    output{i}{k}.time(:,j) = output{i}{k}.time(:,j) - output{i}{k}.time(1,j); % begin time at 0
                else % this handles cases when the first element of the column is a NaN (i.e., no first time)
                    output{i}{k}.time(:,j) = output{i}{k}.time(:,j) - firstTime;
                end
            end
            
            % calculate duration of long windows of nans (>50 ms)
            for j = 1:Ntrials/Nblock
                nans = find(isnan(output{i}{k}.time(:,j)));
                span = diff(nans);
                allOnes = span == 1;
                len = CountOnes(allOnes);
                longDrops{j} = len(len > 3)*(1/frameRate);
            end
            
            % extract whether block was mirrored
            mirror = d.mirror(d.block==k);
            
            for j = 1:length(fields)
                output{i}{k}.(fields{j}) = sineParams.(fields{j})(1:Ntrials/Nblock);
            end
            output{i}{k}.mirror = mirror{1};
            output{i}{k}.trialType = d.trialType(d.block == 1);
            output{i}{k}.frameRate = frameRate;
            output{i}{k}.frameDrops = frameDrops;
            output{i}{k}.longDrops = longDrops;
            output{i}{k}.OS = d.OS{1};
            output{i}{k}.browser = d.browser{1};
%             output{i}.cmConvert = d.cmConvert(1);
        end
    end
end

% compute the number of consecutive ones in an input vector
function len = CountOnes(v)
n   = length(v);
len = zeros(1, ceil(n/2));%, 'uint32');
j   = 0;
k = 1;
while k <= n
  if v(k)
    a = k;
    k = k + 1;
    while k <= n && v(k)
       k = k + 1;
    end
    j = j + 1;
    len(j) = k - a;
  end
  k = k + 1;
end
len = len(1:j);
end