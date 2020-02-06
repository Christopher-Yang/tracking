function data = analyze_data(d, subj_name, block_name, uw)

    disp('Analyzing...');
    rng(3);
    
    % set variables for analysis
    Nsubj = length(subj_name);
    Nblocks = length(block_name);
    Nreps = size(d{1}.(block_name{1}).traj,3);
    names = {'x_x','y_y','x_y','y_x'};
    names2 = {'x_x_all','y_y_all','x_y_all','y_x_all'};
    outputs = {'cursor','Rhand'};
    
    for i = 1:Nsubj
        disp(['   ' subj_name{i}]);
        for j = 1:Nblocks
            output = d{i}.(block_name{j}).traj; % hand, cursor, and target position
            input = d{i}.(block_name{j}).tFile; % frequencies, amplitudes, and phases of target sinusoids
            num = length(input)/6;
            
            freqs = input(1:num*2);
            sorted_freqs = sort(freqs);
            freqs_x = input(1:num);
            freqs_y = input(num+1:num*2);
            amplitudes_x = input(num*2+1:num*3);
            amplitudes_y = input(num*3+1:num*4);
            
            % store hand, cursor, and target position data into trajs or trajs_all
            for k = 0:3
                trajs(:,:,k+1) = [mean(output(:,k*2+1,:),3) mean(output(:,k*2+2,:),3)]';
                for l = 1:Nreps
                    trajs_all(:,:,l,k+1) = ([output(:,k*2+1,l) output(:,k*2+2,l)]);
                end
                trajs(:,:,k+1) = trajs(:,:,k+1) - repmat([mean(trajs(1,:,k+1)) mean(trajs(2,:,k+1))]', [1 size(trajs,2)]); % position averaged across trials
                trajs_all(:,:,:,k+1) = trajs_all(:,:,:,k+1) - repmat(mean(trajs_all(:,:,:,k+1),1), [size(trajs_all,1) 1 1]); % position in all trial
            end
            
            % create data structures to store all data; *_pos stores the
            % position averaged across trials and *_pos_all stores the
            % position for each trial; Rhand stands for right hand
            target = struct('x_pos',trajs(1,:,1)','y_pos',trajs(2,:,1)','x_pos_all',squeeze(trajs_all(:,1,:,1)),'y_pos_all',squeeze(trajs_all(:,2,:,1)));
            Rhand = struct('x_pos',trajs(1,:,3)','y_pos',trajs(2,:,3)','x_pos_all',squeeze(trajs_all(:,1,:,3)),'y_pos_all',squeeze(trajs_all(:,2,:,3)));
            cursor = struct('x_pos',trajs(1,:,4)','y_pos',trajs(2,:,4)','x_pos_all',squeeze(trajs_all(:,1,:,4)),'y_pos_all',squeeze(trajs_all(:,2,:,4)));
            
            % compute mean-squared error between cursor and target for
            % every trial
            MSE = mean((cursor.x_pos_all-target.x_pos_all).^2 + (cursor.y_pos_all-target.y_pos_all).^2,1);
            
            fs = 130.004; % sampling rate for data collection
            x_axis = fs*(0:size(Rhand.x_pos_all,1)/2)/size(Rhand.x_pos_all,1); % frequencies used for plotting amplitude spectra
            time = output(:,11,1)-output(1,11,1); % time at which data was collected
            
            % computes the phasors, gain, and phase between the target and 
            % either the right hand or cursor for every trial
            data{i}.(block_name{j}).phasors.Rhand = fourier2(Rhand,target,freqs_x,freqs_y);
            data{i}.(block_name{j}).phasors.cursor = fourier2(cursor,target,freqs_x,freqs_y);
            
            % average frequency-domain data across trials
            for l = 1:length(outputs)
                for k = 1:length(names)
                    data{i}.(block_name{j}).phasors.(outputs{l}).(names{k}).ratio = mean(data{i}.(block_name{j}).phasors.(outputs{l}).(names2{k}).ratio,2); % average phasors acros trials
                    data{i}.(block_name{j}).phasors.(outputs{l}).(names{k}).gain = abs(data{i}.(block_name{j}).phasors.(outputs{l}).(names{k}).ratio); % compute gain
                    data{i}.(block_name{j}).phasors.(outputs{l}).(names{k}).phase = unwrap(angle(data{i}.(block_name{j}).phasors.(outputs{l}).(names{k}).gain)); % compute phase
                end
            end
            
            % compute amplitude spectra
            [Rhand.x_fft, Rhand.x_fftRaw] = fourier(Rhand.x_pos_all); 
            [Rhand.y_fft, Rhand.y_fftRaw] = fourier(Rhand.y_pos_all);
            [cursor.x_fft, cursor.x_fftRaw] = fourier(cursor.x_pos_all);
            [cursor.y_fft, cursor.y_fftRaw] = fourier(cursor.y_pos_all);
            [target.x_fft, target.x_fftRaw] = fourier(target.x_pos_all);
            [target.y_fft, target.y_fftRaw] = fourier(target.y_pos_all);
            
            % store data
            data{i}.(block_name{j}).cursor = cursor;
            data{i}.(block_name{j}).Rhand = Rhand;
            data{i}.(block_name{j}).target = target;
            data{i}.(block_name{j}).time = time;
            data{i}.(block_name{j}).MSE = MSE;
            for k = 1:length(names)
                all.cursor.(names{k})(j,:,i) = data{i}.(block_name{j}).phasors.cursor.(names{k}).ratio;
                all.Rhand.(names{k})(j,:,i) = data{i}.(block_name{j}).phasors.Rhand.(names{k}).ratio;
            end
        end
    end
    
    % average data acorss subjects
    disp('   averaging...');
    data{Nsubj+1}.freqX = freqs_x;
    data{Nsubj+1}.freqY = freqs_y;
    data{Nsubj+1}.ampX = amplitudes_x;
    data{Nsubj+1}.ampY = amplitudes_y;
    data{Nsubj+1}.x_axis = x_axis;
    
    for k = 1:2
        for i = 1:length(names)
            FT = mean(all.(outputs{k}).(names{i}),3); % take mean of complex ratios
            GAIN = abs(FT);
            PHASE = unwrap(angle(FT),[],2);
            data{Nsubj+1}.(outputs{k}).(names{i}).fft = FT;
            data{Nsubj+1}.(outputs{k}).(names{i}).gain = GAIN;
            data{Nsubj+1}.(outputs{k}).(names{i}).phase = PHASE;
            
            switch names{i}
                case {'x_x','x_y'}
                    data{Nsubj+1}.(outputs{k}).(names{i}).index = data{1}.(block_name{1}).phasors.cursor.x_x_all.index;
                case {'y_y','y_x'}
                    data{Nsubj+1}.(outputs{k}).(names{i}).index = data{1}.(block_name{1}).phasors.cursor.y_y_all.index;
            end
        end
    end
end

% this function simply organizes the data and frequencies and inputs them
% appropriately to fourier()
function out = fourier2(output,input,freqs_x,freqs_y)
    out.x_x_all = fourier(output.x_pos_all,input.x_pos_all,length(freqs_x));
    out.y_y_all = fourier(output.y_pos_all,input.y_pos_all,length(freqs_y));
    out.x_y_all = fourier(output.y_pos_all,input.x_pos_all,length(freqs_x));
    out.y_x_all = fourier(output.x_pos_all,input.y_pos_all,length(freqs_y));
end