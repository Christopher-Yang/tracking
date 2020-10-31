function data = analyze_data(d)
    
    disp('Analyzing...');
    rng(3);
    Nsubj = length(d);
    fields = {'time','cursorX','cursorY','handX','handY','targetX','targetY','cursorX_input','cursorY_input'};
    fields2 = {'tX_freq','tY_freq','tX_amp','tY_amp','tX_phase','tY_phase','cX_freq','cY_freq','cX_amp','cY_amp','cX_phase','cY_phase'};

    for i = 1:Nsubj
        disp(['   Subject ' num2str(i)]);
        Nblock = length(d{i});
        for k = 1:Nblock
            Ntrials = size(d{i}{k}.time,2);
            
            % store sine parmaeters in "input"
            for j = 1:length(fields2)
                input.(fields2{j}) = d{i}{k}.(fields2{j});
            end
            
            % simulate the cursor sinusoidal perturbations
            %         [d{i}.cSimX, d{i}.cSimY] = reconstruct_cursor(input, d{i}.time, d{i}.frameRate);
            
            % linearly interpolate missing data
            %         fields3 = [fields 'cSimX' 'cSimY'];
            for j = 1:length(fields)
                d{i}{k}.(fields{j}) = interpolate_data(d{i}{k}.(fields{j}));
            end
            
            cursorX_raw = d{i}{k}.cursorX + d{i}{k}.cursorX_input;
            cursorY_raw = d{i}{k}.cursorY + d{i}{k}.cursorY_input;
            
            % compute mean-squared error
            MSE = nanmean((cursorX_raw-d{i}{k}.targetX).^2 + (cursorY_raw-d{i}{k}.targetY).^2);
            
            % create data structures to store position data
            trajs.target = struct('x',d{i}{k}.targetX,'y',d{i}{k}.targetY);
            trajs.cursorHand = struct('x',d{i}{k}.handX,'y',d{i}{k}.handY);
            trajs.cursorInput = struct('x',d{i}{k}.cursorX_input,'y',d{i}{k}.cursorY_input);
            
            % calculate frequencies used for FFT
            fs = d{i}{k}.frameRate;
            x_axis = fs*(0:length(trajs.target.x)/2)/length(trajs.target.x);
            
            % perform FFTs and compute complex ratios
            [data{i}{k}.phasors, data{i}{k}.raw_fft, data{i}{k}.processed_fft] = fourier(trajs, input, x_axis, d{i}{k}.trialType);
            
            trajs.rawCursor = struct('x',cursorX_raw,'y',cursorY_raw);
            
            a = data{i}{k}.phasors;
            pAxis = fieldnames(a);
            pAxis = pAxis(1:8);
            Nfreq = length(d{1}{1}.tX_freq{1})*4;
            
            clear Hur Hud 
            for j = 1:4
                if j <= 2
                    for m = 1:2
                        Hur(j,:,m) = reshape(permute([a.(pAxis{j}){2*(m-1)+1}.ratio a.(pAxis{j}){2*(m-1)+2}.ratio a.(pAxis{j}){2*(m-1)+5}.ratio a.(pAxis{j}){2*(m-1)+6}.ratio],[2 1]),[Nfreq 1]);
                        Hud(j,:,m) = reshape(permute([a.(pAxis{j+4}){2*(m-1)+2}.ratio a.(pAxis{j+4}){2*(m-1)+1}.ratio a.(pAxis{j+4}){2*(m-1)+6}.ratio a.(pAxis{j+4}){2*(m-1)+5}.ratio],[2 1]),[Nfreq 1]);
                    end
                else
                    for m = 1:2
                        Hur(j,:,m) = reshape(permute([a.(pAxis{j}){2*(m-1)+5}.ratio a.(pAxis{j}){2*(m-1)+6}.ratio a.(pAxis{j}){2*(m-1)+1}.ratio a.(pAxis{j}){2*(m-1)+2}.ratio],[2 1]),[Nfreq 1]);
                        Hud(j,:,m) = reshape(permute([a.(pAxis{j+4}){2*(m-1)+6}.ratio a.(pAxis{j+4}){2*(m-1)+5}.ratio a.(pAxis{j+4}){2*(m-1)+2}.ratio a.(pAxis{j+4}){2*(m-1)+1}.ratio],[2 1]),[Nfreq 1]);
                    end
                end
            end
            
            Hur = reshape(Hur,[2 2 Nfreq 2]);
            Hud = reshape(Hud,[2 2 Nfreq 2]);
            
            if strcmp(d{i}{k}.mirror,'TRUE') 
                M = [0 1; 1 0];
            else
                M = eye(2);
            end
            
            % testing whether B and F analysis is okay with fixed Hur and
            % Hud
%             interval = pi/6;
%             angle = 0:interval:2*pi-interval;
%             realPart = cos(angle);
%             imaginaryPart = sin(angle);
%             comp = permute(realPart + imaginaryPart*1j,[2 1]);
%             comp = 0.7.*comp;
%             
%             if strcmp(d{i}{k}.mirror,'TRUE') 
%                 Hur(1,2,:) = comp;
%                 Hur(2,1,:) = comp;
%                 Hur = cat(4,Hur,Hur);
% 
%                 Hud(1,2,:) = comp;
%                 Hud(2,1,:) = comp;
%                 Hud = cat(4,Hud,Hud);
%             else
%                 Hur(1,1,:) = comp;
%                 Hur(2,2,:) = comp;
%                 Hur = cat(4,Hur,Hur);
% 
%                 Hud(1,1,:) = comp;
%                 Hud(2,2,:) = comp;
%                 Hud = cat(4,Hud,Hud);
%             end
            
            for m = 1:2
                for j = 1:Nfreq
                    B(:,:,j,m) = -Hud(:,:,j,m)*inv(M*Hud(:,:,j,m) + eye(2));
                    F(:,:,j,m) = Hur(:,:,j,m) + B(:,:,j,m)*(M*Hur(:,:,j,m) - eye(2));
                end
            end
            
            % store data
            data{i}{k}.Hur = Hur;
            data{i}{k}.Hud = Hud;
            data{i}{k}.B = B;
            data{i}{k}.F = F;
            data{i}{k}.pos = trajs;
            data{i}{k}.x_axis = x_axis;
            data{i}{k}.MSE = MSE;
            data{i}{k}.sineParams = input;
            data{i}{k}.trialType = d{i}{k}.trialType;
            data{i}{k}.time = d{i}{k}.time;
            data{i}{k}.frameDrops = d{i}{k}.frameDrops;
            data{i}{k}.longDrops = d{i}{k}.longDrops;
        end
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
                idx = 1;
                while isnan(d(1))
                    num = d(idx+2)-d(idx+1);
                    if ~isnan(num)
                       for i = 0:idx-1
                           d(idx-i) = d(idx-i+1) - num;
                       end
                    end
                    idx = idx + 1;
                end
            elseif isempty(after)
                idx = length(d);
                while isnan(d(end))
                    num = d(idx-2)-d(idx-1);
                    if ~isnan(num)
                        for i = 0:length(d)-idx
                            d(idx+i) = d(idx+i-1) - num;
                        end
                    end
                    idx = idx - 1;
                end
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
