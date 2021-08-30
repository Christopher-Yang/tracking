function graph_lag(data)
    load flipSign
    groups = {'day2','day5','day10'};
%     blocks = {'B1_baseline','B9','B11_habit';
%               'B1_baseline','B24','B26_habit';
%               'B1_baseline','B49','B51_habit'};
    blocks = {'B1_baseline','B9';
              'B1_baseline','B24';
              'B1_baseline','B49'};
    
    f_x = data.day2{1}.B1_baseline.freqX';
    f_y = data.day2{1}.B1_baseline.freqY';
    Nblock = size(blocks,2);
    Nfreq = length(f_x);
    col = [180 180 0
           0 191 255
           255 99 71]./255;

    for i = 1:length(groups)
        Nsubj = length(data.(groups{i}));
        
        for j = 1:Nsubj
            for k = 1:Nblock
                if i == 3 && j == 4 && k == 1
                    phaseX.(groups{i})(:,:,k,j) = [NaN(Nfreq,1) data.(groups{i}){j}.(blocks{i,k}).cursor.phasors.x_x.phase];
                    phaseY.(groups{i})(:,:,k,j) = [NaN(Nfreq,1) data.(groups{i}){j}.(blocks{i,k}).cursor.phasors.y_y.phase];
                else
                    phaseX.(groups{i})(:,:,k,j) = data.(groups{i}){j}.(blocks{i,k}).cursor.phasors.x_x.phase;
                    phaseY.(groups{i})(:,:,k,j) = data.(groups{i}){j}.(blocks{i,k}).cursor.phasors.y_y.phase;
                end
            end
        end
    end
    
    phaseY.day2(3:6,1,2,8) = phaseY.day2(3:6,1,2,8) - 2*pi;
    phaseY.day2(6,3,2,8) = phaseY.day2(6,3,2,8) - 2*pi;
    phaseY.day5(5:6,4,2,5) = phaseY.day5(5:6,4,2,5) - 2*pi;
    phaseX.day5(6,5,2,6) = phaseX.day5(6,5,2,6) - 2*pi;
    phaseX.day5(2:5,5,2,9) = phaseX.day5(2:5,5,2,9) - 2*pi;
    phaseX.day5(5:6,2,2,9) = phaseX.day5(5:6,2,2,9) - 2*pi;
    phaseY.day5(2:6,3,2,9) = phaseY.day5(2:6,3,2,9) - 2*pi;    
    phaseY.day5(5:6,5,2,9) = phaseY.day5(5:6,5,2,9) - 2*pi;    
    phaseY.day10(5:6,[2 4 5],2,1) = phaseY.day10(5:6,[2 4 5],2,1) - 2*pi;    
    phaseY.day10(6,1,2,5) = phaseY.day10(6,1,2,5) - 2*pi;    
    
    for i = 1:length(groups)
        Nsubj = length(data.(groups{i}));
        
        delay.(groups{i}).x = -1000*phaseX.(groups{i})./(repmat(f_x,[1 5 Nblock Nsubj])*2*pi);
        delay.(groups{i}).y = -1000*phaseY.(groups{i})./(repmat(f_y,[1 5 Nblock Nsubj])*2*pi);
        
        dx = diff(mean(delay.(groups{i}).x,2,'omitnan'),[],3);
        dy = diff(mean(delay.(groups{i}).y,2,'omitnan'),[],3);
        
        d = reshape(permute([dx dy],[2 1 4 3]),[12 Nsubj]);
        
        delayMu.(groups{i}) = squeeze(mean(d,2,'omitnan'));        
        delaySE.(groups{i}) = squeeze(std(d,[],2,'omitnan'))./Nsubj;
    end
    
    figure(1); clf; hold on
    for i = 1:length(groups)
        plot([0 0],[1000 1000],'Color',col(i,:),'LineWidth',2)
    end
    plot([0 2],[0 0],'k')
    for i = 1:length(groups)
        s = shadedErrorBar(sort([f_x; f_y]),delayMu.(groups{i}),delaySE.(groups{i}));
        editErrorBar(s,col(i,:),1);        
    end
    set(gca,'TickDir','out')
    xticks(0:2)
    yticks(-250:250:500)
    xlabel('Frequency (Hz)')
    ylabel('Increase in lag (ms)')
    axis([0 2 -400 500])
    legend({'2-day','5-day','10-day'})
end