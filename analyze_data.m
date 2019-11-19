function data = analyze_data(d, block_name, uw, rotate)
    
    disp('Analyzing...');
    rng(3);
    if rotate == 1
        R = rotz(-45);
        R = R(1:2, 1:2);
    else
        R = eye(2);
    end
    Nsubj = length(d);
    Nblocks = length(block_name);
    Nreps = size(d{1}.(block_name{1}).traj,3);
%     names = {'x_x','y_y','x_y','y_x'};
%     names2 = {'x_x_all','y_y_all','x_y_all','y_x_all'};
    outputs = {'cursor'};

    for i = 1:Nsubj
        disp(['   Subject ' num2str(i)]);
        for j = 1:Nblocks
            output = d{i}.(block_name{j}).traj;
            input = d{i}.(block_name{j}).tFile;
            num = length(input)/6;
            
            freqs = input(1:num);
            amplitudes = input(num*2+1:num*3);
            
            for k = 0:3
                trajs(:,:,k+1) = R*([mean(output(:,k*2+1,:),3) mean(output(:,k*2+2,:),3)]');
                for l = 1:Nreps
                    trajs_all(:,:,l,k+1) = (R*([output(:,k*2+1,l) output(:,k*2+2,l)])')';
                end
                trajs(:,:,k+1) = trajs(:,:,k+1) - repmat([mean(trajs(1,:,k+1)) mean(trajs(2,:,k+1))]', [1 size(trajs,2)]);
                trajs_all(:,:,:,k+1) = trajs_all(:,:,:,k+1) - repmat(mean(trajs_all(:,:,:,k+1),1), [size(trajs_all,1) 1 1]);
            end
            
%             create data structures to store all data
            target = struct('x_pos',trajs(1,:,1)','x_pos_all',squeeze(trajs_all(:,1,:,1))); %no time shift
            cursor = struct('x_pos',trajs(1,:,4)','x_pos_all',squeeze(trajs_all(:,1,:,4)));
            
            % compute mean-squared error
            MSE = mean((cursor.x_pos_all-target.x_pos_all).^2);
            
            % compute cross-correlation
%             for k = 1:Nreps
%                 corMat = corrcoef(Lhand.x_pos_all(:,k),Rhand.x_pos_all(:,k));
%                 xCor(k) = corMat(1,2);
%                 corMat = corrcoef(Lhand.y_pos_all(:,k),Rhand.y_pos_all(:,k));
%                 yCor(k) = corMat(1,2);
%             end
            
            fs = 130.004;
            x_axis = fs*(0:length(cursor.x_pos)/2)/length(cursor.x_pos);
            time = output(:,11,1)-output(1,11,1);
            
%             data{i}.(block_name{j}).phasors.Lhand = fourier2(Lhand,target,freqs,freqs_y);
%             data{i}.(block_name{j}).phasors.Rhand = fourier2(Rhand,target,freqs,freqs_y);
            data{i}.(block_name{j}).phasors.cursor = fourier2(cursor,target,freqs);
            
            for l = 1:length(outputs)
                data{i}.(block_name{j}).phasors.cursor.x_x.ratio = mean(data{i}.(block_name{j}).phasors.cursor.x_x_all.ratio,2);
                data{i}.(block_name{j}).phasors.cursor.x_x.gain = abs(data{i}.(block_name{j}).phasors.cursor.x_x.ratio);
                data{i}.(block_name{j}).phasors.cursor.x_x.phase = unwrap(angle(data{i}.(block_name{j}).phasors.cursor.x_x.gain));
            end
            
            cursor.x_fft = fourier(cursor.x_pos_all);   %perform fft and get amplitude data
%             cursor.y_fft = fourier(cursor.y_pos_all);
            target.x_fft = fourier(target.x_pos_all);
%             target.y_fft = fourier(target.y_pos_all);
%             Rhand.x_fft = fourier(Rhand.x_pos_all);
%             Rhand.y_fft = fourier(Rhand.y_pos_all);
%             Lhand.x_fft = fourier(Lhand.x_pos_all);
%             Lhand.y_fft = fourier(Lhand.y_pos_all);
            
            data{i}.(block_name{j}).cursor = cursor;
            data{i}.(block_name{j}).target = target;
%             data{i}.(block_name{j}).Rhand = Rhand;
%             data{i}.(block_name{j}).Lhand = Lhand;
            data{i}.(block_name{j}).time = time;
            data{i}.(block_name{j}).MSE = MSE;
%             data{i}.(block_name{j}).xCor = xCor;
%             data{i}.(block_name{j}).yCor = yCor;
%             for k = 1:length(names)
                all.cursor.x_x(j,:,i) = data{i}.(block_name{j}).phasors.cursor.x_x.ratio;
%                 all.Rhand.(names{k})(j,:,i) = data{i}.(block_name{j}).phasors.Rhand.(names{k}).ratio;
%                 all.Lhand.(names{k})(j,:,i) = data{i}.(block_name{j}).phasors.Lhand.(names{k}).ratio;
%             end
        end
    end
    
    disp('   averaging...');
    n = 1000;
    data{Nsubj+1}.freqX = freqs;
%     data{Nsubj+1}.freqY = freqs_y;
    data{Nsubj+1}.ampX = amplitudes;
%     data{Nsubj+1}.ampY = amplitudes_y;
    data{Nsubj+1}.x_axis = x_axis;
    
    for k = 1
%         for i = 1:length(names)
            FT = mean(all.cursor.x_x,3); % take mean of complex ratios
            AMP = abs(FT);
            PHASE = unwrap(angle(FT),[],2);
            data{Nsubj+1}.cursor.x_x.fft = FT;
            data{Nsubj+1}.cursor.x_x.amplitude = AMP;
            data{Nsubj+1}.cursor.x_x.phase = PHASE;
            
            Nblocks = size(FT,1);
            Nfreqs = size(FT,2);
            % bootstrap data
            boot.x_x.fft = NaN(Nblocks,Nfreqs,n);
            for j = 1:n
                y = datasample(all.cursor.x_x,10,3);
                boot.x_x.fft(:,:,j) = mean(y,3);
            end
            
            data{Nsubj+1}.cursor.x_x.all_amp = abs(all.cursor.x_x);
            
            %         boot.x_x.d = fourier(boot.x_x.fft,0);
            boot.x_x.amplitude = abs(boot.x_x.fft);
            y1 = 20*log10(sort(boot.x_x.amplitude,3));
            boot.x_x.phase = angle(boot.x_x.fft);
            boot.x_x.phase = sort(unwrap(sort(boot.x_x.phase,3),[],2),3);
            
            if uw ~= 0
                switch names{i}
                    case {'x_x'}
                        [y2,a] = error_unwrap(boot.x_x.phase,PHASE,(k-1)*4+1);
                        data{Nsubj+1}.cursor.x_x.phase = a;
                    case {'y_y'}
                        [y2,a] = error_unwrap(boot.x_x.phase,PHASE,(k-1)*4+2);
                        data{Nsubj+1}.cursor.x_x.phase = a;
                    case {'x_y'}
                        [y2,a] = error_unwrap(boot.x_x.phase,PHASE,(k-1)*4+3);
                        data{Nsubj+1}.cursor.x_x.phase = a;
                    case {'y_x'}
                        [y2,a] = error_unwrap(boot.x_x.phase,PHASE,(k-1)*4+4);
                        data{Nsubj+1}.cursor.x_x.phase = a;
                end
            else
                y2 = boot.x_x.phase*(180/pi);
            end
            AMP_ERR = cat(3,y1(:,:,26),y1(:,:,n-25));
            AMP_ERR(:,:,1) = 20*log10(AMP) - AMP_ERR(:,:,1);
            AMP_ERR(:,:,2) = AMP_ERR(:,:,2) - 20*log10(AMP);
            data{Nsubj+1}.cursor.x_x.amp_err = AMP_ERR;
            data{Nsubj+1}.cursor.x_x.amp_err_full = y1;
            
            PHASE_ERR = cat(3,y2(:,:,26),y2(:,:,n-25));
            PHASE_ERR(:,:,1) = PHASE - PHASE_ERR(:,:,1);
            PHASE_ERR(:,:,2) = PHASE_ERR(:,:,2) - (180/pi).*PHASE;
            data{Nsubj+1}.cursor.x_x.phase_err = PHASE_ERR;
            data{Nsubj+1}.cursor.x_x.phase_err_full = y2;
            
%             switch names{i}
%                 case {'x_x','x_y'}
                    data{Nsubj+1}.cursor.x_x.index = data{1}.(block_name{1}).phasors.cursor.x_x_all.index;
%                 case {'y_y','y_x'}
%                     data{Nsubj+1}.cursor.x_x.index = data{1}.(block_name{1}).phasors.cursor.y_y_all.index;
%             end
            optimal_mag = cos(PHASE);
            [x,y] = pol2cart(PHASE,optimal_mag);
            data{Nsubj+1}.cursor.x_x.optimal_fft = complex(x,y);
%         end
    end
end

function out = fourier2(output,input,freqs)
    out.x_x_all = fourier(output.x_pos_all,input.x_pos_all,length(freqs));
end