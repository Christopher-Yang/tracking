function data = analyze_data(d, subj_name, block_name, rotate, uw)
    
    disp('Analyzing...');
    rng(3);
    R = rotz(180);
    R = R(1:2, 1:2);
    Nsubj = length(subj_name);
    Nblocks = length(block_name);
    Nreps = size(d.(subj_name{1}).(block_name{1}).traj,3);
    if rotate == true
        d.most_freq = d.most_freq*2;
    end
    p = NaN(Nblocks,d.most_freq,length(subj_name));
%     r = NaN(length(block_name)*8,d.most_freq,length(subj_name));
    all = struct('x_x',p,'y_y',p,'x_y',p,'y_x',p);
    SRcohere = struct('x_x',p,'y_y',p);
%     SRcohere2 = struct('x_x',r,'y_y',r);
    RRcohere = struct('x_x',p,'y_y',p);
    names = {'x_x','y_y','x_y','y_x'};
    freqs_x = NaN;
    freqs_y = NaN;
    for i = 1:Nsubj
        disp(['   ' subj_name{i}]);
        for j = 1:Nblocks
            output = d.(subj_name{i}).(block_name{j}).traj;
            input = d.(subj_name{i}).(block_name{j}).tFile;
            num = length(input)/6;
            
            if rotate == true
                freqs = input(1:num*2);
                sorted_freqs = sort(freqs);
                freqs_x = sorted_freqs;
                freqs_y = sorted_freqs;
                amplitudes_x = sort(input(num*2+1:num*4),'descend');
                amplitudes_y = sort(input(num*2+2:num*4),'descend');
            else
                freqs = input(1:num*2);
                sorted_freqs = sort(freqs);
                freqs_x = input(1:num);
                freqs_y = input(num+1:num*2);
                amplitudes_x = input(num*2+1:num*3);
                amplitudes_y = input(num*3+1:num*4);
            end
            
            c = [mean(output(1:end-13,7,:),3) mean(output(1:end-13,8,:),3)]';
            t = [mean(output(1:end-13,1,:),3) mean(output(1:end-13,2,:),3)]';
%             c = c - repmat([0.8 0.3]', [1 size(c,2)]);
            t = t - repmat([0.8 0.3]', [1 size(t,2)]);
            c = c - repmat([mean(c(1,:)) mean(c(2,:))]', [1 size(c,2)]); % baseline correction for dark trials
            c_all = [output(1:end-13,7,:) output(1:end-13,8,:)];
            t_all = [output(1:end-13,1,:) output(1:end-13,2,:)];
%             c_all = c_all - repmat([0.8 0.3], [size(c_all,1) 1 size(c_all,3)]);
            t_all = t_all - repmat([0.8 0.3], [size(t_all,1) 1 size(t_all,3)]);
            c_all = c_all - repmat(mean(c_all), [size(c_all,1) 1 1]);  % baseline correction for dark trials
            
            MSE = NaN(1,size(output,3));
            for k = 1:size(output,3)
                MSE(k) = mean((c_all(:,1,k)-t_all(:,1,k)).^2 + (c_all(:,2,k)-t_all(:,2,k)).^2);
            end

%             if rotate == true
%                 c = R*c;
%                 t = R*t;
%                 for k = 1:size(c_all,3)
%                     c_all(:,:,k) = c_all(:,:,k)*R';
%                     t_all(:,:,k) = t_all(:,:,k)*R';
%                 end
%             end
            
%             create data structures to store all data
            cursor = struct('x_pos',c(1,:)','y_pos',c(2,:)'); %no time shift
            target = struct('x_pos',t(1,:)','y_pos',t(2,:)');
            cursor_all = struct('x_pos',squeeze(c_all(:,1,:)),'y_pos',squeeze(c_all(:,2,:)));
            target_all = struct('x_pos',squeeze(t_all(:,1,:)),'y_pos',squeeze(t_all(:,2,:)));
            fs = 130.004;
            x_axis = fs*(0:length(cursor.x_pos)/2)/length(cursor.x_pos);
            time = output(:,11,1)-output(1,11,1);
            
            data.(subj_name{i}).(block_name{j}).x_x = fourier(cursor.x_pos,target.x_pos,length(freqs_x));
            data.(subj_name{i}).(block_name{j}).y_y = fourier(cursor.y_pos,target.y_pos,length(freqs_y));
            data.(subj_name{i}).(block_name{j}).x_y = fourier(cursor.y_pos,target.x_pos,length(freqs_x));
            data.(subj_name{i}).(block_name{j}).y_x = fourier(cursor.x_pos,target.y_pos,length(freqs_y));
            
            for k = 1:length(names)
                data.(subj_name{i}).(block_name{j}).(names{k}).phase = unwrap(data.(subj_name{i}).(block_name{j}).(names{k}).phase)*(180/pi);
            end
            
            data.(subj_name{i}).(block_name{j}).x_x_all = fourier(cursor_all.x_pos, target_all.x_pos, length(freqs_x));
            data.(subj_name{i}).(block_name{j}).y_y_all = fourier(cursor_all.y_pos, target_all.y_pos, length(freqs_y));
            data.(subj_name{i}).(block_name{j}).x_y_all = fourier(cursor_all.y_pos, target_all.x_pos, length(freqs_x));
            data.(subj_name{i}).(block_name{j}).y_x_all = fourier(cursor_all.x_pos, target_all.y_pos, length(freqs_y));
            
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
            
            data.(subj_name{i}).(block_name{j}).x_x.SRcohere = mean(cohxx);
            data.(subj_name{i}).(block_name{j}).y_y.SRcohere = mean(cohyy);
            
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
            
            data.(subj_name{i}).(block_name{j}).x_x.RRcohere = mean(cohx);
            data.(subj_name{i}).(block_name{j}).y_y.RRcohere = mean(cohy);
            
            cursor.x_fft = fourier(cursor.x_pos);   %perform fft and get amplitude data
            cursor.y_fft = fourier(cursor.y_pos);
            target.x_fft = fourier(target.x_pos);
            target.y_fft = fourier(target.y_pos);
            data.(subj_name{i}).(block_name{j}).cursor = cursor;
            data.(subj_name{i}).(block_name{j}).target = target;
            data.(subj_name{i}).(block_name{j}).time = time;
            data.(subj_name{i}).(block_name{j}).MSE = MSE;
            for k = 1:length(names)
                all.(names{k})(j,:,i) = data.(subj_name{i}).(block_name{j}).(names{k}).ratio;
            end
            SRcohere.x_x(j,:,i) = mean(cohxx);
            SRcohere.y_y(j,:,i) = mean(cohyy);
            SRcohere2.x_x(Nreps*(j-1)+1:Nreps*(j-1)+Nreps,:,i) = cohxx;
            SRcohere2.y_y(Nreps*(j-1)+1:Nreps*(j-1)+Nreps,:,i) = cohyy;
            RRcohere.x_x(j,:,i) = data.(subj_name{i}).(block_name{j}).x_x.RRcohere;
            RRcohere.y_y(j,:,i) = data.(subj_name{i}).(block_name{j}).y_y.RRcohere;
        end
    end
    
    disp('   averaging...');
    n = 1000;
    data.avg.x_x.amp = amplitudes_x;
    data.avg.y_y.amp = amplitudes_y;
    data.avg.x_x.x_axis = x_axis;
    
    for i = 1:length(names)
        data.avg.(names{i}).fft = mean(all.(names{i}),3); % take mean of complex ratios
%         data.avg.(names{i}).d = fourier(data.avg.(names{i}).fft,0); % fft averaged complex ratios
        data.avg.(names{i}).amplitude = abs(data.avg.(names{i}).fft);
        data.avg.(names{i}).phase = angle(data.avg.(names{i}).fft);
        data.avg.(names{i}).phase = unwrap(data.avg.(names{i}).phase,[],2); % unwrap phase
        
        Nblocks = size(data.avg.(names{i}).fft,1);
        Nfreqs = size(data.avg.(names{i}).fft,2);
        % bootstrap data
        boot.(names{i}).fft = NaN(Nblocks,Nfreqs,n);
        for j = 1:n
            y = datasample(all.(names{i}),10,3);
            boot.(names{i}).fft(:,:,j) = mean(y,3);
        end
        
%         for j = 1:length(names)
            data.avg.(names{i}).all_amp = abs(all.(names{i}));
            if contains(names{i},'x_y') || contains(names{i},'y_y')
                data.avg.(names{i}).all_amp_neg = -abs(all.(names{i}));
            end
%         end
        
%         boot.(names{i}).d = fourier(boot.(names{i}).fft,0);
        boot.(names{i}).amplitude = abs(boot.(names{i}).fft);
        y1 = 20*log10(sort(boot.(names{i}).amplitude,3));
        boot.(names{i}).phase = angle(boot.(names{i}).fft);
        boot.(names{i}).phase = sort(unwrap(sort(unwrap(sort(unwrap(unwrap(sort(boot.(names{i}).phase,3),[],3),[],2),3),[],3),3),[],3),3);
        
        if uw ~= 0
            switch names{i}
                case {'x_x'}
                    [y2,a] = error_unwrap(boot.(names{i}).phase,data.avg.(names{i}).phase,(uw-1)*4+1);
                    data.avg.(names{i}).phase = a;
                case {'y_y'}
                    [y2,a] = error_unwrap(boot.(names{i}).phase,data.avg.(names{i}).phase,(uw-1)*4+2);
                    data.avg.(names{i}).phase = a;
                case {'x_y'}
                    [y2,a] = error_unwrap(boot.(names{i}).phase,data.avg.(names{i}).phase,(uw-1)*4+3);
                    data.avg.(names{i}).phase = a;
                case {'y_x'}
                    [y2,a] = error_unwrap(boot.(names{i}).phase,data.avg.(names{i}).phase,(uw-1)*4+4);
                    data.avg.(names{i}).phase = a;
            end
        else
            y2 = boot.(names{i}).phase*(180/pi);
        end
        data.avg.(names{i}).amp_err = cat(3,y1(:,:,26),y1(:,:,n-25));
        data.avg.(names{i}).amp_err(:,:,1) = 20*log10(data.avg.(names{i}).amplitude) - data.avg.(names{i}).amp_err(:,:,1);
        data.avg.(names{i}).amp_err(:,:,2) = data.avg.(names{i}).amp_err(:,:,2) - 20*log10(data.avg.(names{i}).amplitude);
        data.avg.(names{i}).amp_err_full = y1;
        
        data.avg.(names{i}).phase_err = cat(3,y2(:,:,26),y2(:,:,n-25));
        data.avg.(names{i}).phase_err(:,:,1) = (180/pi).*data.avg.(names{i}).phase - data.avg.(names{i}).phase_err(:,:,1);
        data.avg.(names{i}).phase_err(:,:,2) = data.avg.(names{i}).phase_err(:,:,2) - (180/pi).*data.avg.(names{i}).phase;
%         for j = 1:numel(data.avg.(names{i}).phase_err)
%             if data.avg.(names{i}).phase_err(j) < 0
%                 data.avg.(names{i}).phase_err(j) = data.avg.(names{i}).phase_err(j) + 360;
%             end
%         end
        data.avg.(names{i}).phase_err_full = y2;
        
        switch names{i}
            case {'x_x','x_y'}
                data.avg.(names{i}).freqs = freqs_x;
                data.avg.(names{i}).index = data.(subj_name{1}).(block_name{1}).x_x.index;
            case {'y_y','y_x'}
                data.avg.(names{i}).freqs = freqs_y;
                data.avg.(names{i}).index = data.(subj_name{1}).(block_name{1}).y_y.index;
        end
        
        if strcmp(names{i},'x_x') || strcmp(names{i},'y_y')
            data.avg.(names{i}).SRcohere_full = SRcohere.(names{i});
            data.avg.(names{i}).SRcohere = mean(SRcohere.(names{i}),3);
            data.avg.(names{i}).SRcohere2 = mean(SRcohere2.(names{i}),3);
            data.avg.(names{i}).RRcohere_full = RRcohere.(names{i});
            data.avg.(names{i}).RRcohere = mean(RRcohere.(names{i}),3);
        end
        optimal_mag = cos(data.avg.(names{i}).phase);
        [x,y] = pol2cart(data.avg.(names{i}).phase,optimal_mag);
        data.avg.(names{i}).optimal_fft = complex(x,y);
    end
    
%     p = NaN(1,6);
%     z = struct('x_x',p,'y_y',p,'x_y',p,'y_x',p);
%     for j = 1:length(names)
%         for i = 1:length(block_name)
% %             z.(names{j})(:,i) = polyfit(data.avg.(names{j}).freqs_x,my_unwrap(unwrap(data.avg.(names{j}).phase(i,:))),1);
%             z.(names{j})(i) = data.avg.(names{j}).freqs'\data.avg.(names{j}).phase(i,:)';
%         end
%         data.avg.(names{j}).phase_slopes = z.(names{j});
%     end
end