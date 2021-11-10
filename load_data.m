% load_data() extracts raw data from data files
% 
%   subj_name: name of subjects
%   block_name: name of blocks
%   folder: file path to data
%   remove: flag to determine how to remove trials with bad data
% 
% The output data structure "output" contains data from each analyzed
% block. The following fields can be found under each block
% 
%   traj: data from 60 seconds of the trial
%   tFile: frequencies, amplitudes and phases of sines used to run the
%          experiment

function output = load_data(subj_name, block_name, folder, remove)
    
    disp('Loading...');
    samplingRate = 130; % sampling rate of the Flock of Birds
    time = 60; % number of seconds to analyze from each file
    Nsamples = round(time*samplingRate); % number of samples to analyze

    for i = 1:length(subj_name) % loop over subjects
        disp(['   Subject ' num2str(i)]);
        
        for j = 1:length(block_name) % loop over blocks
            
            % full file path to data
            path = [folder,subj_name{i},'/',block_name{j}];

            % read frequencies, amplitudes, and phases from tFile used to
            % run experiment
            tFile = dlmread([path,'/tFile.tgt']);
            fnames = dir(path); % file names in data path
            full = NaN(Nsamples, 13, 5); % preallocate matrix
            
            % store data from each file
            for k = 1:length(fnames(not([fnames.isdir])))-1
                data = dlmread([path,'/',fnames(3+(k-1)).name],' ',6,0); % read data
                idx = data(:,9)==1; % find timepoints where the target is moving
                data = data(idx,:); % extract these data points
                data = data(5*samplingRate+1:5*samplingRate+Nsamples,:); % extra data after 5 sec ramp up
                full(:,:,k) = data;
            end
            
            bName = block_name{j}; % name of current block
            
            % handles removal of bad data from 10-day group, subj4, 
            % B1_baseline
            if remove == 1
                if i == 4 && j == 1
                    output{i}.(bName).traj = full(:,:,2:5);
                else
                    output{i}.(bName).traj = full;
                end
                
            % handles removal of bad data from 2-day group, subj1,
            % B2_darkBaseline
            elseif remove == 2
                if i == 1 && j == 2 % 2-day, subj1, B2_darkBaseline
                    output{i}.(bName).traj = full(:,:,2:5);
                else
                    output{i}.(bName).traj = full;
                end
                
            % store movement data
            else
                output{i}.(bName).traj = full;
            end
            
            % store sinusoid's frequencies, amplitudes, and phases
            output{i}.(bName).tFile = tFile;
        end
    end
end