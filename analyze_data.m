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
            Ntrials = size(d{i}.(block_name{j}).traj,3); % number of trials
            output = d{i}.(block_name{j}).traj; % recorded position data
            input = d{i}.(block_name{j}).tFile; % input freqs, amps, and phases
            Nsines = length(input)/6; % number of sines
            
            freqX = input(1:Nsines); % x frequencies
            freqY = input(Nsines+1:Nsines*2); % y frequencies
            ampX = input(Nsines*2+1:Nsines*3); % x amplitudes
            ampY = input(Nsines*3+1:Nsines*4); % y amplitudes
            
            
%             trajs = NaN(size(output,1),2,Ntrials,4); 
%             for k = 0:3
%                 for l = 1:Ntrials
%                     trajs(:,:,l,k+1) = [output(:,k*2+1,l) output(:,k*2+2,l)];
%                 end
%                 
%                 % baseline subtraction
%                 trajs(:,:,:,k+1) = trajs(:,:,:,k+1) - repmat(mean(trajs(:,:,:,k+1),1), [size(trajs,1) 1 1]);
%             end
            
            % create data structures to store position data
            target = struct('x_pos',squeeze(output(:,1,:)),'y_pos',squeeze(output(:,2,:)));
            Lhand = struct('x_pos',squeeze(output(:,3,:)),'y_pos',squeeze(output(:,4,:))); 
            Rhand = struct('x_pos',squeeze(output(:,5,:)),'y_pos',squeeze(output(:,6,:)));
            cursor = struct('x_pos',squeeze(output(:,7,:)),'y_pos',squeeze(output(:,8,:)));
            
            % perform Fourier analysis
            data{i}.(block_name{j}).Lhand.phasors = fourier2(Lhand,target,freqX,freqY);
            data{i}.(block_name{j}).Rhand.phasors = fourier2(Rhand,target,freqX,freqY);
            data{i}.(block_name{j}).cursor.phasors = fourier2(cursor,target,freqX,freqY);
            
%             [Rhand.x_fft, Rhand.x_fftRaw] = fourier(Rhand.x_pos);
%             [Rhand.y_fft, Rhand.y_fftRaw] = fourier(Rhand.y_pos);
%             [Lhand.x_fft, Lhand.x_fftRaw] = fourier(Lhand.x_pos);
%             [Lhand.y_fft, Lhand.y_fftRaw] = fourier(Lhand.y_pos);
%             [cursor.x_fft, cursor.x_fftRaw] = fourier(cursor.x_pos);
%             [cursor.y_fft, cursor.y_fftRaw] = fourier(cursor.y_pos);
%             [target.x_fft, target.x_fftRaw] = fourier(target.x_pos);
%             [target.y_fft, target.y_fftRaw] = fourier(target.y_pos);
            
            % store everything in data
%             fnames = fieldnames(Rhand);
%             for k = 1:length(fnames)
%                 data{i}.(block_name{j}).cursor.(fnames{k}) = cursor.(fnames{k});
%                 data{i}.(block_name{j}).Lhand.(fnames{k}) = Lhand.(fnames{k});
%                 data{i}.(block_name{j}).Rhand.(fnames{k}) = Rhand.(fnames{k});
%                 data{i}.(block_name{j}).target.(fnames{k}) = target.(fnames{k});
%             end
            data{i}.(block_name{j}).freqX = freqX;
            data{i}.(block_name{j}).freqY = freqY;
            data{i}.(block_name{j}).ampX = ampX;
            data{i}.(block_name{j}).ampY = ampY;
        end
    end
end

function out = fourier2(output,input,freqX,freqY)
    out.x_x = fourier(output.x_pos,input.x_pos,length(freqX));
    out.y_y = fourier(output.y_pos,input.y_pos,length(freqY));
    out.x_y = fourier(output.y_pos,input.x_pos,length(freqX));
    out.y_x = fourier(output.x_pos,input.y_pos,length(freqY));
end