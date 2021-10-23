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
            
            freqX = input(1:num);
            freqY = input(num+1:num*2);
            ampX = input(num*2+1:num*3);
            ampY = input(num*3+1:num*4);
            
            trajs = NaN(size(output,1),2,Nreps,4);
            for k = 0:3
%                 trajs(:,:,k+1) = ([mean(output(:,k*2+1,:),3) mean(output(:,k*2+2,:),3)]');
                for l = 1:Nreps
                    trajs(:,:,l,k+1) = [output(:,k*2+1,l) output(:,k*2+2,l)];
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
        end
    end
end

function out = fourier2(output,input,freqX,freqY)
    out.x_x = fourier(output.x_pos,input.x_pos,length(freqX));
    out.y_y = fourier(output.y_pos,input.y_pos,length(freqY));
    out.x_y = fourier(output.y_pos,input.x_pos,length(freqX));
    out.y_x = fourier(output.x_pos,input.y_pos,length(freqY));
end