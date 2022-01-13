% analyze_data() analyzes the data extracted by load_data()
% 
%   d: data structure containing raw data
%   block_name: names of blocks to be analyzed

function data = analyze_data(d, block_name)
    
    disp('Analyzing...');
    rng(3);
    Nsubj = length(d); % number of subjects
    Nblocks = length(block_name); % number of blocks

    for i = 1:Nsubj % loop over subjects
        disp(['   Subject ' num2str(i)]);
        
        for j = 1:Nblocks % loop over blocks
            output = d{i}.(block_name{j}).traj; % recorded position data
            Nsamples = size(output,1); % number of samples per trial
            output(:,1:8,:) = output(:,1:8,:) - repmat(mean(output(:,1:8,:),1), [Nsamples 1 1]); % baseline subtraction
            input = d{i}.(block_name{j}).tFile; % input freqs, amps, and phases
            Nsines = length(input)/6; % number of sines
            
            freqX = input(1:Nsines); % x frequencies
            freqY = input(Nsines+1:Nsines*2); % y frequencies
            ampX = input(Nsines*2+1:Nsines*3); % x amplitudes
            ampY = input(Nsines*3+1:Nsines*4); % y amplitudes
            
            % create data structures to store position data
            target = struct('x_pos',squeeze(output(:,1,:)),'y_pos',squeeze(output(:,2,:)));
            Lhand = struct('x_pos',squeeze(output(:,3,:)),'y_pos',squeeze(output(:,4,:))); 
            Rhand = struct('x_pos',squeeze(output(:,5,:)),'y_pos',squeeze(output(:,6,:)));
            cursor = struct('x_pos',squeeze(output(:,7,:)),'y_pos',squeeze(output(:,8,:)));
            
            % perform Fourier analysis
            data{i}.(block_name{j}).Lhand.phasors = fourier2(Lhand,target,freqX,freqY);
            data{i}.(block_name{j}).Rhand.phasors = fourier2(Rhand,target,freqX,freqY);
            data{i}.(block_name{j}).cursor.phasors = fourier2(cursor,target,freqX,freqY);
            
            MSE = mean((100*(cursor.x_pos-target.x_pos)).^2 + (100*(cursor.y_pos-target.y_pos)).^2,1);
            
            % store everything in data
            fnames = fieldnames(Rhand);
            for k = 1:length(fnames)
                data{i}.(block_name{j}).cursor.(fnames{k}) = cursor.(fnames{k});
                data{i}.(block_name{j}).Lhand.(fnames{k}) = Lhand.(fnames{k});
                data{i}.(block_name{j}).Rhand.(fnames{k}) = Rhand.(fnames{k});
                data{i}.(block_name{j}).target.(fnames{k}) = target.(fnames{k});
            end
            data{i}.(block_name{j}).freqX = freqX;
            data{i}.(block_name{j}).freqY = freqY;
            data{i}.(block_name{j}).ampX = ampX;
            data{i}.(block_name{j}).ampY = ampY;
            data{i}.(block_name{j}).MSE = MSE;
        end
    end
end

function out = fourier2(output,input,freqX,freqY)
    out.x_x = fourier(output.x_pos,input.x_pos,length(freqX));
    out.y_y = fourier(output.y_pos,input.y_pos,length(freqY));
    out.x_y = fourier(output.y_pos,input.x_pos,length(freqX));
    out.y_x = fourier(output.x_pos,input.y_pos,length(freqY));
end