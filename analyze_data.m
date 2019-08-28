function data = analyze_data(d, block_name, uw)
    
    disp('Analyzing...');
    rng(3);
    R = rotz(180);
    R = R(1:2, 1:2);
    Nsubj = length(d);
    Nblocks = length(block_name);
    Nreps = size(d{1}.(block_name{1}).traj,3);
    p = NaN(Nblocks,length(d{1}.(block_name{1}).tFile)/6,length(d));
    SRcohere = struct('x_x',p,'y_y',p);
    RRcohere = struct('x_x',p,'y_y',p);
    names = {'x_x','y_y','x_y','y_x'};
    freqs_x = NaN;
    freqs_y = NaN;
    outputs = {'cursor','Rhand','Lhand'};

    for i = 1:Nsubj
        disp(['   Subject ' num2str(i)]);
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
                trajs_all(:,:,:,k+1) = [output(:,k*2+1,:) output(:,k*2+2,:)];
                if k == 0
                    trajs(:,:,k+1) = trajs(:,:,k+1) - repmat([0.8 0.3]', [1 size(trajs,2)]);
                    trajs_all(:,:,:,k+1) = trajs_all(:,:,:,k+1) - repmat([0.8 0.3], [size(trajs_all,1) 1 size(trajs_all,3)]);
                else
                    trajs(:,:,k+1) = trajs(:,:,k+1) - repmat([mean(trajs(1,:,k+1)) mean(trajs(2,:,k+1))]', [1 size(trajs,2)]);
                    trajs_all(:,:,:,k+1) = trajs_all(:,:,:,k+1) - repmat(mean(trajs_all(:,:,:,k+1),1), [size(trajs_all,1) 1 1]);
                end
            end
            
%             create data structures to store all data
            target = struct('x_pos',trajs(1,:,1)','y_pos',trajs(2,:,1)'); %no time shift
            Lhand = struct('x_pos',trajs(1,:,2)','y_pos',trajs(2,:,2)'); 
            Rhand = struct('x_pos',trajs(1,:,3)','y_pos',trajs(2,:,3)');
            cursor = struct('x_pos',trajs(1,:,4)','y_pos',trajs(2,:,4)');
            
            target_all = struct('x_pos',squeeze(trajs_all(:,1,:,1)),'y_pos',squeeze(trajs_all(:,2,:,1)));
            Lhand_all = struct('x_pos',squeeze(trajs_all(:,1,:,2)),'y_pos',squeeze(trajs_all(:,2,:,2)));
            Rhand_all = struct('x_pos',squeeze(trajs_all(:,1,:,3)),'y_pos',squeeze(trajs_all(:,2,:,3)));
            cursor_all = struct('x_pos',squeeze(trajs_all(:,1,:,4)),'y_pos',squeeze(trajs_all(:,2,:,4)));
            
            MSE = mean((cursor_all.x_pos-target_all.x_pos).^2 + (cursor_all.y_pos-target_all.y_pos).^2);
                
            fs = 130.004;
            x_axis = fs*(0:length(cursor.x_pos)/2)/length(cursor.x_pos);
            time = output(:,11,1)-output(1,11,1);
            
            data{i}.(block_name{j}).phasors.Lhand = fourier2(Lhand,target,Lhand_all,target_all,freqs_x,freqs_y);
            data{i}.(block_name{j}).phasors.Rhand = fourier2(Rhand,target,Rhand_all,target_all,freqs_x,freqs_y);
            data{i}.(block_name{j}).phasors.cursor = fourier2(cursor,target,cursor_all,target_all,freqs_x,freqs_y);
            
            N = size(target.x_pos,1);
            cohxx = NaN(size(target_all.x_pos,2),length(freqs));
            cohyy = NaN(size(target_all.x_pos,2),length(freqs));
            for k = 1:size(cursor_all.x_pos,2)
                coh = mscohere([cursor_all.x_pos(:,k) cursor_all.y_pos(:,k)],[target.x_pos target.y_pos],blackmanharris(round(N/5)),[],sorted_freqs,fs,'mimo')';
%                 coh = mscohere([cursor_all.x_pos cursor_all.y_pos],[target.x_pos target.y_pos],blackmanharris(round(N/5)),[],N,fs,'mimo')';
                cohxx(k,:) = coh(1,:);
                cohyy(k,:) = coh(2,:);
            end

            cohxx = cohxx(:,1:2:end);
            cohyy = cohyy(:,2:2:end);
            
%             cohxx2 = mscohere(target.x_pos,cursor_all.x_pos,blackmanharris(round(N/5)),[],freqs_x,fs);
%             cohyy2 = mscohere(target.y_pos,cursor_all.y_pos,blackmanharris(round(N/5)),[],freqs_y,fs);
            
            data{i}.(block_name{j}).x_x.SRcohere = mean(cohxx);
            data{i}.(block_name{j}).y_y.SRcohere = mean(cohyy);
            
            total = 0;
            for k = 1:size(cursor_all.x_pos,2)-1
                total = total + k;
            end

            FLAG = 1;
            k1 = 1;
            k2 = size(cursor_all.x_pos,2);
            k3 = 1;
            cohx = NaN(total,length(freqs_x));
            cohy = NaN(total,length(freqs_y));
            while FLAG
                cohx(k3,:) = mscohere(cursor_all.x_pos(:,k1),cursor_all.x_pos(:,k2),[],[],freqs_x,fs);
                cohy(k3,:) = mscohere(cursor_all.y_pos(:,k1),cursor_all.y_pos(:,k2),[],[],freqs_y,fs);
                k1 = k1 + 1;
                if k1 == k2
                    k1 = 1;
                    k2 = k2 - 1;
                end
                if k2 == 1
                    FLAG = 0;
                end
                k3 = k3 + 1;
            end
            
            data{i}.(block_name{j}).x_x.RRcohere = mean(cohx);
            data{i}.(block_name{j}).y_y.RRcohere = mean(cohy);
            
            cursor.x_fft = fourier(cursor.x_pos);   %perform fft and get amplitude data
            cursor.y_fft = fourier(cursor.y_pos);
            target.x_fft = fourier(target.x_pos);
            target.y_fft = fourier(target.y_pos);
            Rhand.x_fft = fourier(Rhand.x_pos);
            Rhand.y_fft = fourier(Rhand.y_pos);
            Lhand.x_fft = fourier(Lhand.x_pos);
            Lhand.y_fft = fourier(Lhand.y_pos);
            data{i}.(block_name{j}).cursor = cursor;
            data{i}.(block_name{j}).target = target;
            data{i}.(block_name{j}).Rhand = Rhand;
            data{i}.(block_name{j}).Lhand = Lhand;
            data{i}.(block_name{j}).time = time;
            data{i}.(block_name{j}).MSE = MSE;
            for k = 1:length(names)
                all.cursor.(names{k})(j,:,i) = data{i}.(block_name{j}).phasors.cursor.(names{k}).ratio;
                all.Rhand.(names{k})(j,:,i) = data{i}.(block_name{j}).phasors.Rhand.(names{k}).ratio;
                all.Lhand.(names{k})(j,:,i) = data{i}.(block_name{j}).phasors.Lhand.(names{k}).ratio;
            end
            SRcohere.x_x(j,:,i) = mean(cohxx);
            SRcohere.y_y(j,:,i) = mean(cohyy);
            SRcohere2.x_x(Nreps*(j-1)+1:Nreps*(j-1)+Nreps,:,i) = cohxx;
            SRcohere2.y_y(Nreps*(j-1)+1:Nreps*(j-1)+Nreps,:,i) = cohyy;
            RRcohere.x_x(j,:,i) = data{i}.(block_name{j}).x_x.RRcohere;
            RRcohere.y_y(j,:,i) = data{i}.(block_name{j}).y_y.RRcohere;
        end
    end
    
    disp('   averaging...');
    n = 1000;
    data{Nsubj+1}.amp = amplitudes_x;
    data{Nsubj+1}.amp = amplitudes_y;
    data{Nsubj+1}.x_axis = x_axis;
    
    for k = 1:3
        for i = 1:length(names)
            FT = mean(all.(outputs{k}).(names{i}),3); % take mean of complex ratios
            AMP = abs(FT);
            PHASE = angle(FT);
            PHASE = unwrap(PHASE,[],2); % unwrap phase
            data{Nsubj+1}.(outputs{k}).(names{i}).fft = FT;
            data{Nsubj+1}.(outputs{k}).(names{i}).amplitude = AMP;
            data{Nsubj+1}.(outputs{k}).(names{i}).phase = PHASE;
            
            Nblocks = size(FT,1);
            Nfreqs = size(FT,2);
            % bootstrap data
            boot.(names{i}).fft = NaN(Nblocks,Nfreqs,n);
            for j = 1:n
                y = datasample(all.(outputs{k}).(names{i}),10,3);
                boot.(names{i}).fft(:,:,j) = mean(y,3);
            end
            
            data{Nsubj+1}.(outputs{k}).(names{i}).all_amp = abs(all.(outputs{k}).(names{i}));
            
            %         boot.(names{i}).d = fourier(boot.(names{i}).fft,0);
            boot.(names{i}).amplitude = abs(boot.(names{i}).fft);
            y1 = 20*log10(sort(boot.(names{i}).amplitude,3));
            boot.(names{i}).phase = angle(boot.(names{i}).fft);
            boot.(names{i}).phase = sort(unwrap(sort(boot.(names{i}).phase,3),[],2),3);
            
            if uw ~= 0
                switch names{i}
                    case {'x_x'}
                        [y2,a] = error_unwrap(boot.(names{i}).phase,PHASE,(uw-1)*4+1);
                        data{Nsubj+1}.(outputs{k}).(names{i}).phase = a;
                    case {'y_y'}
                        [y2,a] = error_unwrap(boot.(names{i}).phase,PHASE,(uw-1)*4+2);
                        data{Nsubj+1}.(outputs{k}).(names{i}).phase = a;
                    case {'x_y'}
                        [y2,a] = error_unwrap(boot.(names{i}).phase,PHASE,(uw-1)*4+3);
                        data{Nsubj+1}.(names{i}).phase = a;
                    case {'y_x'}
                        [y2,a] = error_unwrap(boot.(names{i}).phase,PHASE,(uw-1)*4+4);
                        data{Nsubj+1}.(names{i}).phase = a;
                end
            else
                y2 = boot.(names{i}).phase*(180/pi);
            end
            AMP_ERR = cat(3,y1(:,:,26),y1(:,:,n-25));
            AMP_ERR(:,:,1) = 20*log10(AMP) - AMP_ERR(:,:,1);
            AMP_ERR(:,:,2) = AMP_ERR(:,:,2) - 20*log10(AMP);
            data{Nsubj+1}.(outputs{k}).(names{i}).amp_err = AMP_ERR;
            data{Nsubj+1}.(outputs{k}).(names{i}).amp_err_full = y1;
            
            PHASE_ERR = cat(3,y2(:,:,26),y2(:,:,n-25));
            PHASE_ERR(:,:,1) = PHASE - PHASE_ERR(:,:,1);
            PHASE_ERR(:,:,2) = PHASE_ERR(:,:,2) - (180/pi).*PHASE;
            data{Nsubj+1}.(outputs{k}).(names{i}).phase_err = PHASE_ERR;
            data{Nsubj+1}.(outputs{k}).(names{i}).phase_err_full = y2;
            
            switch names{i}
                case {'x_x','x_y'}
                    data{Nsubj+1}.(outputs{k}).(names{i}).freqs = freqs_x;
                    data{Nsubj+1}.(outputs{k}).(names{i}).index = data{1}.(block_name{1}).phasors.cursor.x_x.index;
                case {'y_y','y_x'}
                    data{Nsubj+1}.(outputs{k}).(names{i}).freqs = freqs_y;
                    data{Nsubj+1}.(outputs{k}).(names{i}).index = data{1}.(block_name{1}).phasors.cursor.y_y.index;
            end
            
            if strcmp(names{i},'x_x') || strcmp(names{i},'y_y')
                data{Nsubj+1}.(outputs{k}).(names{i}).SRcohere_full = SRcohere.(names{i});
                data{Nsubj+1}.(outputs{k}).(names{i}).SRcohere = mean(SRcohere.(names{i}),3);
                data{Nsubj+1}.(outputs{k}).(names{i}).SRcohere2 = mean(SRcohere2.(names{i}),3);
                data{Nsubj+1}.(outputs{k}).(names{i}).RRcohere_full = RRcohere.(names{i});
                data{Nsubj+1}.(outputs{k}).(names{i}).RRcohere = mean(RRcohere.(names{i}),3);
            end
            optimal_mag = cos(PHASE);
            [x,y] = pol2cart(PHASE,optimal_mag);
            data{Nsubj+1}.(outputs{k}).(names{i}).optimal_fft = complex(x,y);
        end
    end
end

function out = fourier2(output,input,output_all,input_all,freqs_x,freqs_y)
    out.x_x = fourier(output.x_pos,input.x_pos,length(freqs_x));
    out.y_y = fourier(output.y_pos,input.y_pos,length(freqs_y));
    out.x_y = fourier(output.y_pos,input.x_pos,length(freqs_x));
    out.y_x = fourier(output.x_pos,input.y_pos,length(freqs_y));
    out.x_x_all = fourier(output_all.x_pos,input_all.x_pos,length(freqs_x));
    out.y_y_all = fourier(output_all.y_pos,input_all.y_pos,length(freqs_y));
    out.x_y_all = fourier(output_all.y_pos,input_all.x_pos,length(freqs_x));
    out.y_x_all = fourier(output_all.x_pos,input_all.y_pos,length(freqs_y));
end