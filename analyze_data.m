function data = analyze_data(d)
    
    disp('Analyzing...');
    rng(3);
    Nsubj = length(d);
    fields = {'time','cursorX','cursorY','targetX','targetY','cursorX_input','cursorY_input'};
%     fields = {'time','cursorX','cursorY','targetX','targetY'};
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
            
            cursorX_raw = d{i}{k}.cursorX + d{i}{k}.cursorX;
            cursorY_raw = d{i}{k}.cursorY + d{i}{k}.cursorY;
            
            % compute mean-squared error
            MSE = nanmean((cursorX_raw-d{i}{k}.targetX).^2 + (cursorY_raw-d{i}{k}.targetY).^2);
            
            % create data structures to store position data
            trajs.target = struct('x',d{i}{k}.targetX,'y',d{i}{k}.targetY);
            %         trajs.cursorHand = struct('x',d{i}.cursorX - d{i}.cSimX,'y',d{i}.cursorY - d{i}.cSimY);
            %         trajs.cursorInput = struct('x',d{i}.cSimX,'y',d{i}.cSimY);
            %         trajs.cursorHand = struct('x',d{i}.cursorX - d{i}.cursorX_input,'y',d{i}.cursorY - d{i}.cursorY_input);
            trajs.cursorHand = struct('x',d{i}{k}.cursorX,'y',d{i}{k}.cursorY);
            trajs.cursorInput = struct('x',d{i}{k}.cursorX_input,'y',d{i}{k}.cursorY_input);
            
            % calculate frequencies used for FFT
            fs = d{i}{k}.frameRate;
            x_axis = fs*(0:length(trajs.target.x)/2)/length(trajs.target.x);
            
            % perform FFTs and compute complex ratios
            [data{i}{k}.phasors, data{i}{k}.raw_fft, data{i}{k}.processed_fft] = fourier(trajs, input, x_axis, d{i}{k}.trialType);
            
            trajs.rawCursor = struct('x',cursorX_raw,'y',cursorY_raw);
            
            a = data{i}{k}.phasors;
            b = data{i}{k}.raw_fft;
            pAxis = fieldnames(a);
            pAxis = pAxis(1:8);
            names = {'x_x','x_y','y_x','y_y'};
%             Nfreq = length(d{1}{1}.tX_freq{1})*4;
            Nfreq = length(d{1}{1}.tX_freq{1});
            
            clear Hur Hud Hur2 Hud2
            for j = 1:4
                Hur(j,:) = a.(pAxis{j}){1}.ratio;
                Hud(j,:) = a.(pAxis{j+4}){2}.ratio;
                Hur2(j,:) = reshape(permute([a.(pAxis{j}){3}.ratio a.(pAxis{j}){4}.ratio],[2 1]),[Nfreq 1]);
                Hud2(j,:) = reshape(permute([a.(pAxis{j}){4}.ratio a.(pAxis{j}){3}.ratio],[2 1]),[Nfreq 1]);
                
%                 if j <= 2
%                     for m = 1:2
%                         Hur(j,:,m) = reshape(permute([a.(pAxis{j}){2*(m-1)+1}.ratio a.(pAxis{j}){2*(m-1)+2}.ratio a.(pAxis{j}){2*(m-1)+5}.ratio a.(pAxis{j}){2*(m-1)+6}.ratio],[2 1]),[Nfreq 1]);
%                         Hud(j,:,m) = reshape(permute([a.(pAxis{j+4}){2*(m-1)+2}.ratio a.(pAxis{j+4}){2*(m-1)+1}.ratio a.(pAxis{j+4}){2*(m-1)+6}.ratio a.(pAxis{j+4}){2*(m-1)+5}.ratio],[2 1]),[Nfreq 1]);
%                     end
%                 else
%                     for m = 1:2
%                         Hur(j,:,m) = reshape(permute([a.(pAxis{j}){2*(m-1)+5}.ratio a.(pAxis{j}){2*(m-1)+6}.ratio a.(pAxis{j}){2*(m-1)+1}.ratio a.(pAxis{j}){2*(m-1)+2}.ratio],[2 1]),[Nfreq 1]);
%                         Hud(j,:,m) = reshape(permute([a.(pAxis{j+4}){2*(m-1)+6}.ratio a.(pAxis{j+4}){2*(m-1)+5}.ratio a.(pAxis{j+4}){2*(m-1)+2}.ratio a.(pAxis{j+4}){2*(m-1)+1}.ratio],[2 1]),[Nfreq 1]);
%                     end
%                 end
            end
            
%             Hur = reshape(Hur,[2 2 Nfreq 2]);
%             Hud = reshape(Hud,[2 2 Nfreq 2]);
            
            Hur = reshape(Hur,[2 2 Nfreq]);
            Hud = reshape(Hud,[2 2 Nfreq]);
            Hur2 = reshape(Hur2,[2 2 Nfreq]);
            Hud2 = reshape(Hud2,[2 2 Nfreq]);
            
            
%             if i == 1
%                 if k == 1
%                     M = eye(2);
%                 else
                    M = [0 1; 1 0];
%                 end
%             if i == 1
%                 M = [0 1; 1 0];
%             elseif i == 2
%                 M = eye(2);
%             end

            for m = 1
                for j = 1:Nfreq
%                 B(:,:,j) = inv(Hud(:,:,j)*M + eye(2))*-Hud(:,:,j);
%                 F(:,:,j) = Hur(:,:,j) + Hur(:,:,j)*M*B(:,:,j) - B(:,:,j);
%                 B2(:,:,j) = inv(Hud2(:,:,j)*M + eye(2))*-Hud2(:,:,j);
%                 F2(:,:,j) = Hur2(:,:,j) + Hur2(:,:,j)*M*B2(:,:,j) - B2(:,:,j);
%                 
                    B(:,:,j,m) = -Hud(:,:,j,m)*M*inv(Hud(:,:,j,m) + eye(2));
                    F(:,:,j,m) = Hur(:,:,j,m) + B(:,:,j,m)*(M*Hur(:,:,j,m) - eye(2));
                    B2(:,:,j,m) = -Hud2(:,:,j,m)*M*inv(Hud2(:,:,j,m) + eye(2));
                    F2(:,:,j,m) = Hur2(:,:,j,m) + B2(:,:,j,m)*(M*Hur2(:,:,j,m) - eye(2));
                end
            end
            
%             for j = 1:length(pAxis)/2
%                 Hur = a.(pAxis{j}){1}.ratio;
%                 Hud = a.(pAxis{j+4}){2}.ratio;
%                 
%                 Hur2 = reshape(permute([a.(pAxis{j}){3}.ratio a.(pAxis{j}){4}.ratio],[2 1]),[Nfreq 1]);
%                 Hud2 = reshape(permute([a.(pAxis{j+4}){4}.ratio a.(pAxis{j+4}){3}.ratio],[2 1]),[Nfreq 1]);
%                 
%                 B.(names{j}) = -Hud./(1+Hud);
%                 F.(names{j}) = (Hur + Hud)./(1 + Hud);
%                 
%                 B2.(names{j}) = -Hud2./(1+Hud2);
%                 F2.(names{j}) = (Hur2 + Hud2)./(1 + Hud2);
%             end
            
            % store data
            data{i}{k}.Hur = Hur;
            data{i}{k}.Hud = Hud;
            data{i}{k}.B = B;
            data{i}{k}.F = F;
            data{i}{k}.B2 = B2;
            data{i}{k}.F2 = F2;
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

% simulate cursor sinusoidal perturbations from sine parameters
function [outputX, outputY] = reconstruct_cursor(input, time, frameRate)
    Ntrials = size(time,2);
    outputX = NaN(size(time));
    outputY = NaN(size(time));
    
    for i = 1:Ntrials
        t = 0:1/frameRate:45-(1/frameRate);
        window = 5*frameRate:45*frameRate-1;
        x = repmat(input.cX_amp{i}',[1 length(t)]).*cos(2*pi*input.cX_freq{i}'*t + repmat(input.cX_phase{i}',[1 length(t)]));
        y = repmat(input.cY_amp{i}',[1 length(t)]).*cos(2*pi*input.cY_freq{i}'*t + repmat(input.cY_phase{i}',[1 length(t)]));
        x = sum(x(:,window),1);
        y = sum(y(:,window),1);

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
