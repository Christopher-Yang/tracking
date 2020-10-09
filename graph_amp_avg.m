function graph_amp_avg(data,gblocks,subj,trialType)

col = [30 144 255
    0 0 139
    255 182 193]./255;
lw = 1;

Ntrials = length(gblocks);
Nsubj = length(data);
fields = {'tX_freq','tY_freq','cX_freq','cY_freq'};
fields2 = {'tX_amp','tY_amp','cX_amp','cY_amp'};
fields3 = {'tX_idx','tY_idx','cX_idx','cY_idx'};

for i = 1:Ntrials
    for j = 1:length(fields)
        input.(fields{j}) = data{1}{gblocks(i)}.sineParams.(fields{j}){trialType};
        input.(fields2{j}) = data{1}{gblocks(i)}.sineParams.(fields2{j}){trialType};
        input.(fields3{j}) = data{1}{gblocks(i)}.phasors.index.(fields{j}){trialType};
    end
    
    xAxis = data{1}{gblocks(i)}.x_axis;
%     trialType = data{1}{gblocks(i)}.trialType;
    
%     Nsteps = length(data{1}.(block_name{1}).(output).x_pos);
    
    g = NaN(length(xAxis),Nsubj);
    inX = g;
    inY = g;
    inX2 = g;
    inY2 = g;
    outX = g;
    outY = g;
    
    for j = 1:Nsubj
        dat = data{j}{gblocks(i)}.processed_fft;
        if trialType == 1
            inX(:,j) = dat.target.xFFT(:,trialType);
            inY(:,j) = dat.target.yFFT(:,trialType);
        elseif trialType == 2
            inX(:,j) = dat.cursorInput.xFFT(:,trialType);
            inY(:,j) = dat.cursorInput.yFFT(:,trialType);
        else
            inX(:,j) = dat.target.xFFT(:,trialType);
            inY(:,j) = dat.target.yFFT(:,trialType);
            inX2(:,j) = dat.cursorInput.xFFT(:,trialType);
            inY2(:,j) = dat.cursorInput.yFFT(:,gblocks(i));
        end
        outX(:,j) = dat.cursorHand.xFFT(:,trialType);
        outY(:,j) = dat.cursorHand.yFFT(:,trialType);
    end
    
    if isempty(subj)
        inX_mean = mean(inX,2);
        inY_mean = mean(inY,2);
        inX2_mean = mean(inX2,2);
        inY2_mean = mean(inY2,2);
        outX_mean = mean(outX,2);
        outY_mean = mean(outY,2);
    else
        inX_mean = inX(:,subj);
        inY_mean = inY(:,subj);
        inX2_mean = inX2(:,subj);
        inY2_mean = inY2(:,subj);
        outX_mean = outX(:,subj);
        outY_mean = outY(:,subj);
    end
    
    if trialType == 1
        morePlots = 0;
        name = 'Target sines';
        ix = input.tX_idx;
        iy = input.tY_idx;
        fx = input.tX_freq;
        fy = input.tY_freq;
    elseif trialType == 2
        morePlots = 0;
        name = 'Cursor sines';
        ix = input.cX_idx;
        iy = input.cY_idx;
        fx = input.cX_freq;
        fy = input.cY_freq;
    else
        morePlots = 1;
        name = 'Target sines';
        ix = input.tX_idx;
        iy = input.tY_idx;
        fx = input.tX_freq;
        fy = input.tY_freq;
        
        name2 = 'Cursor sines';
        ix2 = input.cX_idx;
        iy2 = input.cY_idx;
        fx2 = input.cX_freq;
        fy2 = input.cY_freq;
    end
    
    figure(i); clf
    subplot(2,2,1); hold on
    plot(xAxis,inX_mean,'k')
    plot(xAxis,outX_mean,'Color',col(1,:),'LineWidth',lw)
    plot(fx,outX_mean(ix),'ko-','LineWidth',lw,'MarkerFaceColor',[0 0 0])
    plot(fy,outX_mean(iy),'ko-','LineWidth',lw)
    set(gca,'box','off','TickDir','out')
    axis([0 2 0 2])
    ylabel('X Amplitude (m)')
    title(name)
    legend({'Input','Output','X freqs','Y freqs'})
    legend boxoff
    
    subplot(2,2,3); hold on
    plot(xAxis,inY_mean,'k')
    plot(xAxis,outY_mean,'LineWidth',lw,'Color',col(1,:))
    plot(fy,outY_mean(iy),'ko-','LineWidth',lw)
    plot(fx,outY_mean(ix),'ko-','LineWidth',lw,'MarkerFaceColor',[0 0 0])
    set(gca,'box','off','TickDir','out')
    axis([0 2 0 2])
    ylabel('Y Amplitude (m)')
    xlabel('Frequency (Hz)')
    
    if morePlots
        subplot(2,2,2); hold on
        plot(xAxis,inX2_mean,'k')
        plot(xAxis,outX_mean,'Color',col(1,:),'LineWidth',lw)
        plot(fx2,outX_mean(ix2),'ko-','LineWidth',lw,'MarkerFaceColor',[0 0 0])
        plot(fy2,outX_mean(iy2),'ko-','LineWidth',lw)
        set(gca,'box','off','TickDir','out')
        axis([0 2 0 2])
        title(name2)
        
        subplot(2,2,4); hold on
        plot(xAxis,inY2_mean,'k')
        plot(xAxis,outY_mean,'LineWidth',lw,'Color',col(1,:))
        plot(fy2,outY_mean(iy2),'ko-','LineWidth',lw)
        plot(fx2,outY_mean(ix2),'ko-','LineWidth',lw,'MarkerFaceColor',[0 0 0])
        set(gca,'box','off','TickDir','out')
        axis([0 2 0 2])
        xlabel('Frequency (Hz)')
    end
end
end