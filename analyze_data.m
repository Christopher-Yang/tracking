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

function data = analyze_data(d, subj_name, block_name)

disp('Analyzing...');
rng(3);

% set variables for analysis
Nsubj = length(subj_name);
Nblocks = length(block_name);
Ntrials = size(d{1}.(block_name{1}).traj,3);
Nfreq = length(d{1}.(block_name{1}).tFile)/3;
% names = {'x_x','y_y','x_y','y_x'};
% outputs = {'cursor','Rhand'};

for i = 1:Nsubj % iterate over subjects
    disp(['   ' subj_name{i}]);
    for j = 1:Nblocks % iterate over blocks
        % hand, cursor, and target position
        output = d{i}.(block_name{j}).traj;
        
        % frequencies, amplitudes, and phases of target sinusoids
        input = d{i}.(block_name{j}).tFile;
        
        trialType = d{i}.(block_name{j}).trialType;
        bimanual_mode = d{i}.(block_name{j}).bimanual_mode;
        rotation = d{i}.(block_name{j}).rotation;
                
        % Store hand, cursor, and target position from every trial into 
        % trajs.
        for k = 0:3
            for l = 1:Ntrials
                trajs(:,:,l,k+1) = ([output(:,k*2+1,l) ...
                    output(:,k*2+2,l)]);
            end
            trajs(:,:,:,k+1) = trajs(:,:,:,k+1) - ...
                repmat(mean(trajs(:,:,:,k+1),1), ...
                [size(trajs,1) 1 1]); % apply baseline correction
        end
        
        % reorganize trajs into separate variables for the right hand 
        % (Rhand), cursor, and target
        target = struct('x',squeeze(trajs(:,1,:,1)), ...
            'y',squeeze(trajs(:,2,:,1)));
        Rhand = struct('x',squeeze(trajs(:,1,:,3)), ...
            'y',squeeze(trajs(:,2,:,3)));
        cursor = struct('x',squeeze(trajs(:,1,:,4)), ...
            'y',squeeze(trajs(:,2,:,4)));
        
        rotMat = rotz(-rotation);
        rotMat = rotMat(1:2,1:2);
        
        if bimanual_mode == 0
            mirMat = eye(2);
        elseif bimanual_mode == 3
            mirMat = [0 1; 1 0];
        end
        
        Rhand_rot = permute(trajs(:,:,:,3),[2 1 3]);
        for k = 1:size(Rhand_rot,3)
            Rhand_rot(:,:,k) = mirMat*rotMat*Rhand_rot(:,:,k);
        end
        
        cursorHand = struct('x',squeeze(Rhand_rot(1,:,:)), ...
            'y',squeeze(Rhand_rot(2,:,:)));
        cursorInput = struct('x',cursor.x - cursorHand.x, ...
            'y',cursor.y - cursorHand.y);
        
        trajs2.target = target;
        trajs2.cursorHand = cursorHand;
        trajs2.cursorInput = cursorInput;
        
        % compute mean-squared error between cursor and target for
        % every trial
        MSE = mean((100*(cursor.x-target.x)).^2 + ...
            (100*(cursor.y-target.y)).^2,1);
                
        fs = 130; % sampling rate for data collection
        
        % frequencies used for plotting amplitude spectra
        x_axis = fs*(0:size(Rhand.x,1)/2)/size(Rhand.x,1);
        
        % time at which data was collected
        time = output(:,11,1)-output(1,11,1);
        
        % perform FFTs
        [data{i}.(block_name{j}).phasors, data{i}.(block_name{j}).raw_fft, ...
            data{i}.(block_name{j}).processed_fft] = fourier(trajs2, trialType);

        a = data{i}.(block_name{j}).phasors;
        clear Hur Hud
        for k = 1:Ntrials/4
            m = (k-1)*4;
            
            Hur(1,1,:,k) = reshape(permute([a.xTarg_x{1+m}.ratio a.xTarg_x{2+m}.ratio a.xTarg_x{3+m}.ratio a.xTarg_x{4+m}.ratio], [2 1]), [12 1]);
            Hur(2,1,:,k) = reshape(permute([a.xTarg_y{1+m}.ratio a.xTarg_y{2+m}.ratio a.xTarg_y{3+m}.ratio a.xTarg_y{4+m}.ratio], [2 1]), [12 1]);
            Hur(1,2,:,k) = reshape(permute([a.yTarg_x{4+m}.ratio a.yTarg_x{4+m}.ratio a.yTarg_x{1+m}.ratio a.yTarg_x{2+m}.ratio], [2 1]), [12 1]);
            Hur(2,2,:,k) = reshape(permute([a.yTarg_y{4+m}.ratio a.yTarg_y{4+m}.ratio a.yTarg_y{1+m}.ratio a.yTarg_y{2+m}.ratio], [2 1]), [12 1]);
            
            Hud(1,1,:,k) = reshape(permute([a.xCurs_x{3+m}.ratio a.xCurs_x{4+m}.ratio a.xCurs_x{1+m}.ratio a.xCurs_x{2+m}.ratio], [2 1]), [12 1]);
            Hud(2,1,:,k) = reshape(permute([a.xCurs_y{3+m}.ratio a.xCurs_y{4+m}.ratio a.xCurs_y{1+m}.ratio a.xCurs_y{2+m}.ratio], [2 1]), [12 1]);
            Hud(1,2,:,k) = reshape(permute([a.yCurs_x{2+m}.ratio a.yCurs_x{3+m}.ratio a.yCurs_x{4+m}.ratio a.yCurs_x{1+m}.ratio], [2 1]), [12 1]);
            Hud(2,2,:,k) = reshape(permute([a.yCurs_y{2+m}.ratio a.yCurs_y{3+m}.ratio a.yCurs_y{4+m}.ratio a.yCurs_y{1+m}.ratio], [2 1]), [12 1]);
        end

        M = mirMat*rotMat;
        
        for k = 1:Ntrials/4
            for n = 1:Nfreq
                B(:,:,n,k) = -Hud(:,:,n,k)*inv(M*Hud(:,:,n,k) + eye(2));
                F(:,:,n,k) = Hur(:,:,n,k) + B(:,:,n,k)*(M*Hur(:,:,n,k) - eye(2));
            end
        end
        
        fnames = fieldnames(Rhand);
        
        % store everything in data
        for k = 1:length(fnames)
            data{i}.(block_name{j}).cursor.(fnames{k}) = cursor ...
                .(fnames{k});
            data{i}.(block_name{j}).Rhand.(fnames{k}) = Rhand.(fnames{k});
            data{i}.(block_name{j}).target.(fnames{k}) = target ...
                .(fnames{k});
            data{i}.(block_name{j}).cursorHand.(fnames{k}) = cursorHand ...
                .(fnames{k});
            data{i}.(block_name{j}).cursorInput.(fnames{k}) = cursorInput ...
                .(fnames{k});
        end
        data{i}.(block_name{j}).Hur = Hur;
        data{i}.(block_name{j}).Hud = Hud;
        data{i}.(block_name{j}).F = F;
        data{i}.(block_name{j}).B = B;        
        data{i}.(block_name{j}).time = time;
        data{i}.(block_name{j}).MSE = MSE;
        data{i}.(block_name{j}).x_axis = x_axis;
    end
end
end