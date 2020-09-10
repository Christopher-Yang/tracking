function graph_amp_avg(data,gblocks,output)
    
    col = [30 144 255
           0 0 139
           255 182 193]./255;
    lw = 1;
    
%     H = gobjects(length(gblocks),1);
%     for i = 1:length(H)
%         H(i) = figure;
%     end
    Ntrials = length(gblocks);
    Nsubj = length(data);
    fields = {'tX_freq','tY_freq','cX_freq','cY_freq'};
    fields2 = {'tX_amp','tY_amp','cX_amp','cY_amp'};
    fields3 = {'tX_idx','tY_idx','cX_idx','cY_idx'};
    
    for i = 1:Ntrials
        for j = 1:length(fields)
            input.(fields{j}) = data{1}.sineParams.(fields{j}){gblocks(i)};
            input.(fields2{j}) = data{1}.sineParams.(fields2{j}){gblocks(i)};
            input.(fields3{j}) = data{1}.phasors.index.(fields{j}){gblocks(i)};
        end
        
        xAxis = data{1}.x_axis;
        trialType = data{1}.trialType;
        
%     Nsteps = length(data{1}.(block_name{1}).(output).x_pos);
    
        g = NaN(length(xAxis),Nsubj);
        inX = g;
        inY = g;
        inX2 = g;
        inY2 = g;
        outX = g;
        outY = g;
%         g2 = NaN(length(xAxis),1000,length(gblocks));
%         amps.x = g1;
%         amps.y = g1;
%         a = amps;
    
%     names1 = {'x','y'};
%     names2 = {'bootx','booty'};
%     names3 = {'errx','erry'};
%     names4 = {'x_all','y_all'};
    
        for j = 1:Nsubj
            dat = data{j}.processed_fft;
            if trialType(gblocks(i)) == 1
                inX(:,j) = dat.target.xFFT(:,gblocks(i));
                inY(:,j) = dat.target.yFFT(:,gblocks(i));
            elseif trialType(gblocks(i)) == 2
                inX(:,j) = dat.cursorInput.xFFT(:,gblocks(i));
                inY(:,j) = dat.cursorInput.yFFT(:,gblocks(i));
            else
                inX(:,j) = dat.target.xFFT(:,gblocks(i));
                inY(:,j) = dat.target.yFFT(:,gblocks(i));
                inX2(:,j) = dat.cursorInput.xFFT(:,gblocks(i));
                inY2(:,j) = dat.cursorInput.yFFT(:,gblocks(i));
            end
            outX(:,j) = dat.cursorHand.xFFT(:,gblocks(i));
            outY(:,j) = dat.cursorHand.yFFT(:,gblocks(i));
        end
        
        
%         for j = 1:length(gblocks)
%             dat = data{i}.(block_name{gblocks(j)});
%             pX = dat.(output).x_fft;
%             pY = dat.(output).y_fft;
%             a.x_all(:,j,i) = squeeze(abs(mean(pX,2)));
%             a.y_all(:,j,i) = squeeze(abs(mean(pY,2)));
%         end
    
%     a.x = mean(a.x_all,3);
%     a.y = mean(a.y_all,3);
%     
%     a.x_all = permute(a.x_all,[1 3 2]);
%     a.y_all = permute(a.y_all,[1 3 2]);
    
    % bootstrap data
%     for j = 1:length(names1)
%         for i = 1:1000
%             a.(names2{j})(:,i,:) = mean(datasample(a.(names4{j}),Nsubj,2),2);
%         end
%         a.(names2{j}) = sort(a.(names2{j}),2);
%         a.(names3{j})(:,:,1) = a.(names1{j}) - squeeze(a.(names2{j})(:,26,:));
%         a.(names3{j})(:,:,2) = squeeze(a.(names2{j})(:,975,:)) - a.(names1{j});
%         a.(names3{j}) = permute(a.(names3{j}),[3 1 2]);
%     end
    inX_mean = mean(inX,2);
    inY_mean = mean(inY,2);
    inX2_mean = mean(inX2,2);
    inY2_mean = mean(inY2,2);
    outX_mean = mean(outX,2);
    outY_mean = mean(outY,2);
    
%     for i = 1:length(gblocks)
%         figure(H(i))
    if trialType(gblocks(i)) == 1
        num = 1;
        ix = input.tX_idx;
        iy = input.tY_idx;
    elseif trialType(gblocks(i)) == 2
        num = 1;
        ix = input.cX_idx;
        iy = input.cY_idx;
    else
        num = 2;
        ix = input.tX_idx;
        iy = input.tY_idx;
        ix2 = input.cX_idx;
        iy2 = input.cY_idx;
    end 
    
        figure(i); clf
        subplot(2,num,1); hold on
%         plot(tX_freq,ampX,'ko','LineWidth',lw,'MarkerFaceColor',[0 0 0])
%         plot(tY_freq,ampY,'ko','LineWidth',lw)
        % plot mean with error bars; doesn't work if amplitude spectra are averaged in the complex domain
%         s = shadedErrorBar(xAxis, a.x(:,i),a.errx(:,:,i));
%         editErrorBar(s,col(3,:),0.25);
        plot(xAxis,inX_mean,'k')
        plot(xAxis,outX_mean,'Color',col(1,:),'LineWidth',lw)
        plot(input.tX_freq,outX_mean(input.tX_idx),'ko-','LineWidth',lw,'MarkerFaceColor',[0 0 0])
        plot(input.tY_freq,outX_mean(input.tY_idx),'ko-','LineWidth',lw)
        set(gca,'Xtick',0:1:2,'Ytick',0:0.01:0.03,'box','off','LineWidth',1,'TickDir','out','XMinorTick','on')
        ax = gca;
        ax.XAxis.MinorTickValues = 0:0.25:2.25;
        axis([0 2.3 0 0.025])
        ylabel('X Amplitude (m)')
%         title(graph_name{gblocks(i)})
        legend({'X target','Y target','Hand'})
        legend boxoff
        
        subplot(2,num,2); hold on
%         plot(tX_freq, ampX,'ko','LineWidth',lw,'MarkerFaceColor',[0 0 0])
%         plot(tY_freq, ampY,'ko','LineWidth',lw)
%         s = shadedErrorBar(xAxis, a.y(:,i),a.erry(:,:,i),'lineProps','-b');
%         editErrorBar(s,col(3,:),0.25);
        plot(xAxis,inY_mean,'k')
        plot(xAxis,outY_mean,'LineWidth',lw,'Color',col(1,:))
        plot(input.tY_freq,outY_mean(input.tY_idx),'ko-','LineWidth',lw)
        plot(input.tX_freq,outY_mean(input.tX_idx),'ko-','LineWidth',lw,'MarkerFaceColor',[0 0 0])
        set(gca,'Xtick',0:1:2,'Ytick',0:0.01:0.03,'box','off','LineWidth',1,'TickDir','out','XMinorTick','on')
        ax = gca;
        ax.XAxis.MinorTickValues = 0:0.25:2.25;
        axis([0 2.3 0 0.025])
        ylabel('Y Amplitude (m)')
        xlabel('Frequency (Hz)')
    end
end