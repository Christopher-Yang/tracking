function data = analyze_data(d, subj_name, block_name, uw)
    
    disp('Analyzing...');
    rng(3);
    Nsubj = length(subj_name);
    Nblocks = length(block_name);
    Nreps = size(d{1}.(block_name{1}).traj,3);
    names = {'x_x','y_y','x_y','y_x'};
    names2 = {'x_x_all','y_y_all','x_y_all','y_x_all'};
    outputs = {'cursor','Rhand'};
    for i = 1:Nsubj
        disp(['   ' subj_name{i}]);
        for j = 1:Nblocks
            output = d{i}.(block_name{j}).traj;
            input = d{i}.(block_name{j}).tFile;
            num = length(input)/6;
            
            freqs = input(1:num*2);
            sorted_freqs = sort(freqs);
            freqs_x = input(1:num);
            freqs_y = input(num+1:num*2);
            amplitudes_x = input(num*2+1:num*3);
            amplitudes_y = input(num*3+1:num*4);
            
            for k = 0:3
                trajs(:,:,k+1) = [mean(output(:,k*2+1,:),3) mean(output(:,k*2+2,:),3)]';
                for l = 1:Nreps
                    trajs_all(:,:,l,k+1) = ([output(:,k*2+1,l) output(:,k*2+2,l)]);
                end
                trajs(:,:,k+1) = trajs(:,:,k+1) - repmat([mean(trajs(1,:,k+1)) mean(trajs(2,:,k+1))]', [1 size(trajs,2)]);
                trajs_all(:,:,:,k+1) = trajs_all(:,:,:,k+1) - repmat(mean(trajs_all(:,:,:,k+1),1), [size(trajs_all,1) 1 1]);
            end
            
            %             create data structures to store all data
            target = struct('x_pos',trajs(1,:,1)','y_pos',trajs(2,:,1)','x_pos_all',squeeze(trajs_all(:,1,:,1)),'y_pos_all',squeeze(trajs_all(:,2,:,1))); %no time shift
            Rhand = struct('x_pos',trajs(1,:,3)','y_pos',trajs(2,:,3)','x_pos_all',squeeze(trajs_all(:,1,:,3)),'y_pos_all',squeeze(trajs_all(:,2,:,3)));
            cursor = struct('x_pos',trajs(1,:,4)','y_pos',trajs(2,:,4)','x_pos_all',squeeze(trajs_all(:,1,:,4)),'y_pos_all',squeeze(trajs_all(:,2,:,4)));
            
            MSE = mean((cursor.x_pos_all-target.x_pos_all).^2 + (cursor.y_pos_all-target.y_pos_all).^2,1);
            
            fs = 130.004;
            x_axis = fs*(0:length(Rhand.x_pos)/2)/length(Rhand.x_pos);
            time = output(:,11,1)-output(1,11,1);

            data{i}.(block_name{j}).phasors.Rhand = fourier2(Rhand,target,freqs_x,freqs_y);
            data{i}.(block_name{j}).phasors.cursor = fourier2(cursor,target,freqs_x,freqs_y);
            
            for l = 1:length(outputs)
                for k = 1:length(names)
                    data{i}.(block_name{j}).phasors.(outputs{l}).(names{k}).ratio = mean(data{i}.(block_name{j}).phasors.(outputs{l}).(names2{k}).ratio,2);
                    data{i}.(block_name{j}).phasors.(outputs{l}).(names{k}).gain = abs(data{i}.(block_name{j}).phasors.(outputs{l}).(names{k}).ratio);
                    data{i}.(block_name{j}).phasors.(outputs{l}).(names{k}).phase = unwrap(angle(data{i}.(block_name{j}).phasors.(outputs{l}).(names{k}).gain));
                end
            end
            
            Rhand.x_fft = fourier(Rhand.x_pos_all);   %perform fft and get amplitude data
            Rhand.y_fft = fourier(Rhand.y_pos_all);
            cursor.x_fft = fourier(cursor.x_pos_all);
            cursor.y_fft = fourier(cursor.y_pos_all);
            target.x_fft = fourier(target.x_pos_all);
            target.y_fft = fourier(target.y_pos_all);
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
    
    disp('   averaging...');
    n = 1000;
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
            
            Nblocks = size(FT,1);
            Nfreqs = size(FT,2);
            % bootstrap data
            boot.(names{i}).fft = NaN(Nblocks,Nfreqs,n);
            for j = 1:n
                y = datasample(all.(outputs{k}).(names{i}),10,3);
                boot.(names{i}).fft(:,:,j) = mean(y,3);
            end
            
            data{Nsubj+1}.(outputs{k}).(names{i}).all_gain = abs(all.(outputs{k}).(names{i}));
            data{Nsubj+1}.(outputs{k}).(names{i}).all_phase = angle(all.(outputs{k}).(names{i}));
            
            %         boot.(names{i}).d = fourier(boot.(names{i}).fft,0);
            boot.(names{i}).gain = abs(boot.(names{i}).fft);
            y1 = 20*log10(sort(boot.(names{i}).gain,3));
            boot.(names{i}).phase = angle(boot.(names{i}).fft);
            boot.(names{i}).phase = sort(unwrap(sort(boot.(names{i}).phase,3),[],2),3);
            
            if uw ~= 0
                switch names{i}
                    case {'x_x'}
                        [y2,a] = error_unwrap(boot.(names{i}).phase,PHASE,(k-1)*4+1);
                        data{Nsubj+1}.(outputs{k}).(names{i}).phase = a;
                    case {'y_y'}
                        [y2,a] = error_unwrap(boot.(names{i}).phase,PHASE,(k-1)*4+2);
                        data{Nsubj+1}.(outputs{k}).(names{i}).phase = a;
                    case {'x_y'}
                        [y2,a] = error_unwrap(boot.(names{i}).phase,PHASE,(k-1)*4+3);
                        data{Nsubj+1}.(outputs{k}).(names{i}).phase = a;
                    case {'y_x'}
                        [y2,a] = error_unwrap(boot.(names{i}).phase,PHASE,(k-1)*4+4);
                        data{Nsubj+1}.(outputs{k}).(names{i}).phase = a;
                end
            else
                y2 = boot.(names{i}).phase*(180/pi);
            end
            GAIN_ERR = cat(3,y1(:,:,26),y1(:,:,n-25));
            GAIN_ERR(:,:,1) = 20*log10(GAIN) - GAIN_ERR(:,:,1);
            GAIN_ERR(:,:,2) = GAIN_ERR(:,:,2) - 20*log10(GAIN);
            data{Nsubj+1}.(outputs{k}).(names{i}).gain_err = GAIN_ERR;
            data{Nsubj+1}.(outputs{k}).(names{i}).gain_err_full = y1;
            
            PHASE_ERR = cat(3,y2(:,:,26),y2(:,:,n-25));
            PHASE_ERR(:,:,1) = PHASE - PHASE_ERR(:,:,1);
            PHASE_ERR(:,:,2) = PHASE_ERR(:,:,2) - (180/pi).*PHASE;
            data{Nsubj+1}.(outputs{k}).(names{i}).phase_err = PHASE_ERR;
            data{Nsubj+1}.(outputs{k}).(names{i}).phase_err_full = y2;
            
            switch names{i}
                case {'x_x','x_y'}
                    data{Nsubj+1}.(outputs{k}).(names{i}).index = data{1}.(block_name{1}).phasors.cursor.x_x_all.index;
                case {'y_y','y_x'}
                    data{Nsubj+1}.(outputs{k}).(names{i}).index = data{1}.(block_name{1}).phasors.cursor.y_y_all.index;
            end
            optimal_mag = cos(PHASE);
            [x,y] = pol2cart(PHASE,optimal_mag);
            data{Nsubj+1}.(outputs{k}).(names{i}).optimal_fft = complex(x,y);
        end
    end
end

function out = fourier2(output,input,freqs_x,freqs_y)
    out.x_x_all = fourier(output.x_pos_all,input.x_pos_all,length(freqs_x));
    out.y_y_all = fourier(output.y_pos_all,input.y_pos_all,length(freqs_y));
    out.x_y_all = fourier(output.y_pos_all,input.x_pos_all,length(freqs_x));
    out.y_x_all = fourier(output.x_pos_all,input.y_pos_all,length(freqs_y));
end