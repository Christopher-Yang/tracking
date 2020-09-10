function data = analyze_data(d)
    
    disp('Analyzing...');
    rng(3);
    Nsubj = length(d);
    fields = {'time','cursorX','cursorY','targetX','targetY'};
    fields2 = {'tX_freq','tY_freq','tX_amp','tY_amp','tX_phase','tY_phase','cX_freq','cY_freq','cX_amp','cY_amp','cX_phase','cY_phase'};

    for i = 1:Nsubj
        disp(['   Subject ' num2str(i)]);
        Ntrials = size(d{i}.time,2);
        
        % store sine parmaeters in "input"
        for j = 1:length(fields2)
            input.(fields2{j}) = d{i}.(fields2{j});
        end
        
        % simulate the cursor sinusoidal perturbations 
        [d{i}.cSimX, d{i}.cSimY] = reconstruct_cursor(input, d{i}.time, d{i}.frameRate, d{i}.cmConvert);
        
        % linearly interpolate missing data
        fields3 = [fields 'cSimX' 'cSimY'];
        for j = 1:length(fields3)
            d{i}.(fields3{j}) = interpolate_data(d{i}.(fields3{j}));
        end            
        
        % compute mean-squared error
        MSE = nanmean((d{i}.cursorX-d{i}.targetX).^2 + (d{i}.cursorY-d{i}.targetY).^2);
        
        % create data structures to store position data
        trajs.target = struct('x',d{i}.targetX,'y',d{i}.targetY);
        trajs.cursorHand = struct('x',d{i}.cursorX - d{i}.cSimX,'y',d{i}.cursorY - d{i}.cSimY);
        trajs.cursorInput = struct('x',d{i}.cSimX,'y',d{i}.cSimY);
        
        % calculate frequencies used for FFT
        fs = d{i}.frameRate;
        x_axis = fs*(0:length(trajs.target.x)/2)/length(trajs.target.x);
        
        % perform FFTs and compute complex ratios
        [data{i}.phasors, data{i}.raw_fft, data{i}.processed_fft] = fourier(trajs, input, x_axis, d{i}.trialType);
        
        trajs.rawCursor = struct('x',d{i}.cursorX,'y',d{i}.cursorY);
        
        % store data
        data{i}.pos = trajs;
        data{i}.x_axis = x_axis;
        data{i}.MSE = MSE;
        data{i}.sineParams = input;
        data{i}.trialType = d{i}.trialType;
        data{i}.time = d{i}.time;
        
    end
end

% simulate cursor sinusoidal perturbations from sine parameters
function [outputX, outputY] = reconstruct_cursor(input, time, frameRate, cmConvert)
    Ntrials = size(time,2);
    outputX = NaN(size(time));
    outputY = NaN(size(time));
    
    for i = 1:Ntrials
        t = 0:1/frameRate:40-(1/frameRate);
        x = repmat(input.cX_amp{i}',[1 length(t)]).*cos(2*pi*input.cX_freq{i}'*t + repmat(input.cX_phase{i}',[1 length(t)]));
        y = repmat(input.cY_amp{i}',[1 length(t)]).*cos(2*pi*input.cY_freq{i}'*t + repmat(input.cY_phase{i}',[1 length(t)]));
        x = cmConvert .* sum(x,1);
        y = cmConvert .* sum(y,1);

        idx = find(isnan(time(:,i)));
        x(idx) = NaN;
        y(idx) = NaN;
        
        outputX(:,i) = x;
        outputY(:,i) = y;
    end
end

% interpolates missing data (note that this method may interpolate negative
% times instead of 0 at the start of the trial)
function output = interpolate_data(data)
    Ntrials = size(data,2);
    for k = 1:Ntrials
        d = data(:,k);
        while sum(isnan(d))
            nans = find(isnan(d),1);
            numbers = find(~isnan(d));
            before = numbers(numbers < nans);
            after = numbers(numbers > nans);
            if isempty(before)
                d(1) = d(2) - (d(3)-d(2));
            elseif isempty(after)
                d(end) = d(end-1) - (d(end-2)-d(end-1));
            else
                before = before(end);
                after = after(1);
                d(before:after) = linspace(d(before),d(after),after-before+1);
            end
        end
        data(:,k) = d;
    end
    output = data;
end
