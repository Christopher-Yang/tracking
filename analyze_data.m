function data = analyze_data(d)
    
    disp('Analyzing...');
    rng(3);
    Nsubj = length(d);

    for i = 1:Nsubj
        disp(['   Subject ' num2str(i)]);
        Nreps = size(d{i}.time,2);
%         for j = 1:Nblocks
%             output = d{i}.traj;
            input.cX = d{i}.cX_input;
            input.cY = d{i}.cY_input;
            input.tX = d{i}.tX_input;
            input.tY = d{i}.tY_input;
%             num = length(input)/6;
            
%             freqs = input(1:num*2);
%             sorted_freqs = sort(freqs);
%             freqX = input(1:num);
%             freqY = input(num+1:num*2);
%             ampX = input(num*2+1:num*3);
%             ampY = input(num*3+1:num*4);
            
%             for k = 0:1
%                 for l = 1:Nreps
%                     trajs(:,:,l,k+1) = (([output(:,k*2+2,l) output(:,k*2+3,l)])')';
%                 end
%                 trajs(:,:,:,k+1) = trajs(:,:,:,k+1) - repmat(mean(trajs(:,:,:,k+1),1), [size(trajs,1) 1 1]);
%             end
            
%             create data structures to store all data
            target = struct('x_pos',d{i}.targetX,'y_pos',d{i}.targetY);
            cursor = struct('x_pos',d{i}.cursorX,'y_pos',d{i}.cursorY);
            
            % compute mean-squared error
            MSE = mean((cursor.x_pos-target.x_pos).^2 + (cursor.y_pos-target.y_pos).^2);
            
            fs = 60;
            x_axis = fs*(0:length(cursor.x_pos)/2)/length(cursor.x_pos);
            time = d{i}.time;
            
            data{i}.cursor.phasors = fourier2(cursor,target,input,d{i}.cursorPerturb);
            
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
            fnames = fieldnames(cursor);
            for k = 1:length(fnames)
                data{i}.cursor.(fnames{k}) = cursor.(fnames{k});
                data{i}.target.(fnames{k}) = target.(fnames{k});
            end
            data{i}.time = time;
            data{i}.MSE = MSE;
            data{i}.freqX = freqX;
            data{i}.freqY = freqY;
            data{i}.ampX = ampX;
            data{i}.ampY = ampY;
%             for k = 1:length(names)
%                 all.cursor.(names{k})(j,:,i) = data{i}.(block_name{j}).phasors.cursor.(names{k}).ratio;
%                 all.Rhand.(names{k})(j,:,i) = data{i}.(block_name{j}).phasors.Rhand.(names{k}).ratio;
%                 all.Lhand.(names{k})(j,:,i) = data{i}.(block_name{j}).phasors.Lhand.(names{k}).ratio;
%             end
%         end
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

function out = fourier2(output,input,input,cursorPerturb)
    out.x_x = fourier(output.x_pos,input.x_pos,length(freqX));
    out.y_y = fourier(output.y_pos,input.y_pos,length(freqY));
    out.x_y = fourier(output.y_pos,input.x_pos,length(freqX));
    out.y_x = fourier(output.x_pos,input.y_pos,length(freqY));
end