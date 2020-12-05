function data = analyze_data(d, subj_name, block_name)
% Process raw tracking data, "d", for the given subjects in "subj_name" 
% and blocks in "block_name." The output data structure, "data," is 
% organized as "data1.[group]{[subject number]}.[block]." In the block 
% field, both time- and frequency-domain information for the right hand, 
% cursor, and target movements are stored in the fields "Rhand," "cursor," 
% and "target," respectively. "Rhand" and "cursor" contain the phasors
% between each of these inputs and target movement. For example, 
% phasors.x_y is the phasors between x-target movement and y-hand/cursor 
% movement. The fields containing single matrices have the following 
% contents: 
% 
%   time - timestamp of data collection
%   MSE - mean-squared error for every trial
%   freqX/freqY - x- and y-target frequencies
%   ampX/ampY - x- and y-target amplitudes
%   x_axis - frequencies for plotting amplitude spectra
%   x_pos/y_pos - x- and y-positions for every trial
%   x_fft/y_fft - x and y single-sided amplitude spectra
%   x_fftRaw/y_fftRaw - raw double-sided DFT 
%   ratio - complex ratio (phasor) between input and output
%   gain/phase - gain/phase of the phasor
%   index - indices in x_axis which are the target frequencies

disp('Analyzing...');
rng(3);

% set variables for analysis
Nsubj = length(subj_name);
Nblocks = length(block_name);
Nreps = size(d{1}.(block_name{1}).traj,3);
names = {'x_x','y_y','x_y','y_x'};
outputs = {'cursor','Rhand'};

for i = 1:Nsubj % iterate over subjects
    disp(['   ' subj_name{i}]);
    for j = 1:Nblocks % iterate over blocks
        % hand, cursor, and target position
        output = d{i}.(block_name{j}).traj;
        
        % frequencies, amplitudes, and phases of target sinusoids
        input = d{i}.(block_name{j}).tFile;
        num = length(input)/6; % break up input into 6 sections
        
        % store frequencies and amplitudes in variables
        freqs = input(1:num*2);
        sorted_freqs = sort(freqs);
        freqX = input(1:num);
        freqY = input(num+1:num*2);
        ampX = input(num*2+1:num*3);
        ampY = input(num*3+1:num*4);
        
        % Store hand, cursor, and target position from every trial into 
        % trajs.
        for k = 0:3
            for l = 1:Nreps
                trajs(:,:,l,k+1) = ([output(:,k*2+1,l) ...
                    output(:,k*2+2,l)]);
            end
            trajs(:,:,:,k+1) = trajs(:,:,:,k+1) - ...
                repmat(mean(trajs(:,:,:,k+1),1), ...
                [size(trajs,1) 1 1]); % apply baseline correction
        end
        
        % reorganize trajs into separate variables for the right hand 
        % (Rhand), cursor, and target
        target = struct('x_pos',squeeze(trajs(:,1,:,1)), ...
            'y_pos',squeeze(trajs(:,2,:,1)));
        Rhand = struct('x_pos',squeeze(trajs(:,1,:,3)), ...
            'y_pos',squeeze(trajs(:,2,:,3)));
        cursor = struct('x_pos',squeeze(trajs(:,1,:,4)), ...
            'y_pos',squeeze(trajs(:,2,:,4)));
        
        % compute mean-squared error between cursor and target for
        % every trial
        MSE = mean((100*(cursor.x_pos-target.x_pos)).^2 + ...
            (100*(cursor.y_pos-target.y_pos)).^2,1);
                
        fs = 130.004; % sampling rate for data collection
        
        % frequencies used for plotting amplitude spectra
        x_axis = fs*(0:size(Rhand.x_pos,1)/2)/size(Rhand.x_pos,1);
        
        % time at which data was collected
        time = output(:,11,1)-output(1,11,1);
        
        % computes the phasors, gain, and phase between the target and
        % either the right hand or cursor for every trial
        data{i}.(block_name{j}).Rhand.phasors = fourier2(Rhand, ...
            target,freqX,freqY);
        data{i}.(block_name{j}).cursor.phasors = fourier2(cursor, ...
            target,freqX,freqY);
        
        % compute amplitude spectra
        [Rhand.x_fft, Rhand.x_fftRaw] = fourier(Rhand.x_pos);
        [Rhand.y_fft, Rhand.y_fftRaw] = fourier(Rhand.y_pos);
        [cursor.x_fft, cursor.x_fftRaw] = fourier(cursor.x_pos);
        [cursor.y_fft, cursor.y_fftRaw] = fourier(cursor.y_pos);
        [target.x_fft, target.x_fftRaw] = fourier(target.x_pos);
        [target.y_fft, target.y_fftRaw] = fourier(target.y_pos);
        
        fnames = fieldnames(Rhand);
        
        % store everything in data
        for k = 1:length(fnames)
            data{i}.(block_name{j}).cursor.(fnames{k}) = cursor ...
                .(fnames{k});
            data{i}.(block_name{j}).Rhand.(fnames{k}) = Rhand.(fnames{k});
            data{i}.(block_name{j}).target.(fnames{k}) = target ...
                .(fnames{k});
        end
        data{i}.(block_name{j}).time = time;
        data{i}.(block_name{j}).MSE = MSE;
        data{i}.(block_name{j}).freqX = freqX;
        data{i}.(block_name{j}).freqY = freqY;
        data{i}.(block_name{j}).ampX = ampX;
        data{i}.(block_name{j}).ampY = ampY;
        data{i}.(block_name{j}).x_axis = x_axis;
    end
end

% this function simply organizes the data and frequencies and inputs them
% appropriately to fourier()
function out = fourier2(output,input,freqX,freqY)
    out.x_x = fourier(output.x_pos,input.x_pos,length(freqX));
    out.y_y = fourier(output.y_pos,input.y_pos,length(freqY));
    out.x_y = fourier(output.y_pos,input.x_pos,length(freqX));
    out.y_x = fourier(output.x_pos,input.y_pos,length(freqY));
end
end