function data = analyze_data(d, block_name)
    
    disp('Analyzing...');
    rng(3);
    Nsubj = length(d);
    Nblocks = length(block_name);

    for i = 1:Nsubj
        disp(['   Subject ' num2str(i)]);
        for j = 1:Nblocks
            Nreps = size(d{i}.(block_name{j}).traj,3);
            output = d{i}.(block_name{j}).traj;
            input = d{i}.(block_name{j}).tFile;
            num = length(input)/6;
            
            freqs = input(1:num*2);
            sorted_freqs = sort(freqs);
            freqX = input(1:num);
            freqY = input(num+1:num*2);
            ampX = input(num*2+1:num*3);
            ampY = input(num*3+1:num*4);
            
            for k = 0:3
%                 trajs(:,:,k+1) = ([mean(output(:,k*2+1,:),3) mean(output(:,k*2+2,:),3)]');
                for l = 1:Nreps
                    trajs(:,:,l,k+1) = (([output(:,k*2+1,l) output(:,k*2+2,l)])')';
                end
%                 trajs(:,:,k+1) = trajs(:,:,k+1) - repmat([mean(trajs(1,:,k+1)) mean(trajs(2,:,k+1))]', [1 size(trajs,2)]);
                trajs(:,:,:,k+1) = trajs(:,:,:,k+1) - repmat(mean(trajs(:,:,:,k+1),1), [size(trajs,1) 1 1]);
            end
            
%             create data structures to store all data
            target = struct('x_pos',squeeze(trajs(:,1,:,1)),'y_pos',squeeze(trajs(:,2,:,1)));
            Lhand = struct('x_pos',squeeze(trajs(:,1,:,2)),'y_pos',squeeze(trajs(:,2,:,2))); 
            Rhand = struct('x_pos',squeeze(trajs(:,1,:,3)),'y_pos',squeeze(trajs(:,2,:,3)));
            cursor = struct('x_pos',squeeze(trajs(:,1,:,4)),'y_pos',squeeze(trajs(:,2,:,4)));
            
            % compute mean-squared error
            MSE = mean((cursor.x_pos-target.x_pos).^2 + (cursor.y_pos-target.y_pos).^2);
            
            % compute cross-correlation
            for k = 1:Nreps
                corMat = corrcoef(Lhand.x_pos(:,k),Rhand.x_pos(:,k));
                xCor(k) = corMat(1,2);
                corMat = corrcoef(Lhand.y_pos(:,k),Rhand.y_pos(:,k));
                yCor(k) = corMat(1,2);
            end
            
            fs = 130.004;
            x_axis = fs*(0:length(cursor.x_pos)/2)/length(cursor.x_pos);
            time = output(:,11,1)-output(1,11,1);
            
            data{i}.(block_name{j}).Lhand.phasors = fourier2(Lhand,target,freqX,freqY);
            data{i}.(block_name{j}).Rhand.phasors = fourier2(Rhand,target,freqX,freqY);
            data{i}.(block_name{j}).cursor.phasors = fourier2(cursor,target,freqX,freqY);
            
            [Rhand.x_fft, Rhand.x_fftRaw] = fourier(Rhand.x_pos);
            [Rhand.y_fft, Rhand.y_fftRaw] = fourier(Rhand.y_pos);
            [Lhand.x_fft, Lhand.x_fftRaw] = fourier(Lhand.x_pos);
            [Lhand.y_fft, Lhand.y_fftRaw] = fourier(Lhand.y_pos);
            [cursor.x_fft, cursor.x_fftRaw] = fourier(cursor.x_pos);
            [cursor.y_fft, cursor.y_fftRaw] = fourier(cursor.y_pos);
            [target.x_fft, target.x_fftRaw] = fourier(target.x_pos);
            [target.y_fft, target.y_fftRaw] = fourier(target.y_pos);
            
            
%             for l = 1:length(outputs)
%                 for k = 1:length(names)
%                     data{i}.(block_name{j}).phasors.(outputs{l}).(names{k}).ratio = mean(data{i}.(block_name{j}).phasors.(outputs{l}).(names2{k}).ratio,2);
%                     data{i}.(block_name{j}).phasors.(outputs{l}).(names{k}).gain = abs(data{i}.(block_name{j}).phasors.(outputs{l}).(names{k}).ratio);
%                     data{i}.(block_name{j}).phasors.(outputs{l}).(names{k}).phase = unwrap(angle(data{i}.(block_name{j}).phasors.(outputs{l}).(names{k}).gain));
%                 end
%             end
            
            % store everything in data
            fnames = fieldnames(Rhand);
            for k = 1:length(fnames)
                data{i}.(block_name{j}).cursor.(fnames{k}) = cursor.(fnames{k});
                data{i}.(block_name{j}).Lhand.(fnames{k}) = Lhand.(fnames{k});
                data{i}.(block_name{j}).Rhand.(fnames{k}) = Rhand.(fnames{k});
                data{i}.(block_name{j}).target.(fnames{k}) = target.(fnames{k});
            end
            data{i}.(block_name{j}).time = time;
            data{i}.(block_name{j}).MSE = MSE;
            data{i}.(block_name{j}).xCor = xCor;
            data{i}.(block_name{j}).yCor = yCor;
            data{i}.(block_name{j}).freqX = freqX;
            data{i}.(block_name{j}).freqY = freqY;
            data{i}.(block_name{j}).ampX = ampX;
            data{i}.(block_name{j}).ampY = ampY;
%             for k = 1:length(names)
%                 all.cursor.(names{k})(j,:,i) = data{i}.(block_name{j}).phasors.cursor.(names{k}).ratio;
%                 all.Rhand.(names{k})(j,:,i) = data{i}.(block_name{j}).phasors.Rhand.(names{k}).ratio;
%                 all.Lhand.(names{k})(j,:,i) = data{i}.(block_name{j}).phasors.Lhand.(names{k}).ratio;
%             end
        end
    end
    
%     disp('   averaging...');
%     n = 1000;
%     data{Nsubj+1}.freqX = freqX;
%     data{Nsubj+1}.freqY = freqY;
%     data{Nsubj+1}.ampX = amplitudes_x;
%     data{Nsubj+1}.ampY = amplitudes_y;
%     data{Nsubj+1}.x_axis = x_axis;
    
%     for k = 1:3
%         for i = 1:length(names)
%             FT = mean(all.(outputs{k}).(names{i}),3); % take mean of complex ratios
%             AMP = abs(FT);
%             PHASE = unwrap(angle(FT),[],2);
%             data{Nsubj+1}.(outputs{k}).(names{i}).fft = FT;
%             data{Nsubj+1}.(outputs{k}).(names{i}).amplitude = AMP;
%             data{Nsubj+1}.(outputs{k}).(names{i}).phase = PHASE;
%             
%             Nblocks = size(FT,1);
%             Nfreqs = size(FT,2);
%             % bootstrap data
%             boot.(names{i}).fft = NaN(Nblocks,Nfreqs,n);
%             for j = 1:n
%                 y = datasample(all.(outputs{k}).(names{i}),10,3);
%                 boot.(names{i}).fft(:,:,j) = mean(y,3);
%             end
%             
%             data{Nsubj+1}.(outputs{k}).(names{i}).all_amp = abs(all.(outputs{k}).(names{i}));
%             
%             %         boot.(names{i}).d = fourier(boot.(names{i}).fft,0);
%             boot.(names{i}).amplitude = abs(boot.(names{i}).fft);
%             y1 = 20*log10(sort(boot.(names{i}).amplitude,3));
%             boot.(names{i}).phase = angle(boot.(names{i}).fft);
%             boot.(names{i}).phase = sort(unwrap(sort(boot.(names{i}).phase,3),[],2),3);
%             
%             if uw ~= 0
%                 switch names{i}
%                     case {'x_x'}
%                         [y2,a] = error_unwrap(boot.(names{i}).phase,PHASE,(k-1)*4+1);
%                         data{Nsubj+1}.(outputs{k}).(names{i}).phase = a;
%                     case {'y_y'}
%                         [y2,a] = error_unwrap(boot.(names{i}).phase,PHASE,(k-1)*4+2);
%                         data{Nsubj+1}.(outputs{k}).(names{i}).phase = a;
%                     case {'x_y'}
%                         [y2,a] = error_unwrap(boot.(names{i}).phase,PHASE,(k-1)*4+3);
%                         data{Nsubj+1}.(outputs{k}).(names{i}).phase = a;
%                     case {'y_x'}
%                         [y2,a] = error_unwrap(boot.(names{i}).phase,PHASE,(k-1)*4+4);
%                         data{Nsubj+1}.(outputs{k}).(names{i}).phase = a;
%                 end
%             else
%                 y2 = boot.(names{i}).phase*(180/pi);
%             end
%             AMP_ERR = cat(3,y1(:,:,26),y1(:,:,n-25));
%             AMP_ERR(:,:,1) = 20*log10(AMP) - AMP_ERR(:,:,1);
%             AMP_ERR(:,:,2) = AMP_ERR(:,:,2) - 20*log10(AMP);
%             data{Nsubj+1}.(outputs{k}).(names{i}).amp_err = AMP_ERR;
%             data{Nsubj+1}.(outputs{k}).(names{i}).amp_err_full = y1;
%             
%             PHASE_ERR = cat(3,y2(:,:,26),y2(:,:,n-25));
%             PHASE_ERR(:,:,1) = PHASE - PHASE_ERR(:,:,1);
%             PHASE_ERR(:,:,2) = PHASE_ERR(:,:,2) - (180/pi).*PHASE;
%             data{Nsubj+1}.(outputs{k}).(names{i}).phase_err = PHASE_ERR;
%             data{Nsubj+1}.(outputs{k}).(names{i}).phase_err_full = y2;
%             
%             switch names{i}
%                 case {'x_x','x_y'}
%                     data{Nsubj+1}.(outputs{k}).(names{i}).index = data{1}.(block_name{1}).phasors.cursor.x_x_all.index;
%                 case {'y_y','y_x'}
%                     data{Nsubj+1}.(outputs{k}).(names{i}).index = data{1}.(block_name{1}).phasors.cursor.y_y_all.index;
%             end
%             optimal_mag = cos(PHASE);
%             [x,y] = pol2cart(PHASE,optimal_mag);
%             data{Nsubj+1}.(outputs{k}).(names{i}).optimal_fft = complex(x,y);
%         end
%     end
end

function out = fourier2(output,input,freqX,freqY)
    out.x_x = fourier(output.x_pos,input.x_pos,length(freqX));
    out.y_y = fourier(output.y_pos,input.y_pos,length(freqY));
    out.x_y = fourier(output.y_pos,input.x_pos,length(freqX));
    out.y_x = fourier(output.x_pos,input.y_pos,length(freqY));
end