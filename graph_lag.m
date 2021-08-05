function graph_lag(data)
    load flipSign
    groups = {'day2','day5','day10'};
    names = {'2-day','5-day','10-day'};
    blocks = {'B1_baseline','B9','B11_habit';
              'B1_baseline','B24','B26_habit';
              'B1_baseline','B49','B51_habit'};
    
    f_x = data.day2{1}.B1_baseline.freqX';
    f_y = data.day2{1}.B1_baseline.freqY';
    Nblock = size(blocks,2);
    Nfreq = length(f_x);
    col = copper;
    n = round(64.*((1:Nfreq-1)/(Nfreq-1)));
    col = col([1 n],:);

    for i = 1:length(groups)
        Nsubj = length(data.(groups{i}));
        
        for j = 1:Nsubj
            data.(groups{i}){j}.(blocks{i,3}).cursor.phasors.x_x.phase(flipSign.(groups{i})(:,:,j)) = ...
                data.(groups{i}){j}.(blocks{i,3}).cursor.phasors.x_x.phase(flipSign.(groups{i})(:,:,j)) - pi;
            
            data.(groups{i}){j}.(blocks{i,3}).cursor.phasors.x_x.phase = ...
                unwrap(data.(groups{i}){j}.(blocks{i,3}).cursor.phasors.x_x.phase);
        end
        
%         for j = 1:Nsubj
%             for k = 1:3
%                 phaseX.(groups{i})(:,k,j) = mean(data.(groups{i}){j}.(blocks{i,k}).cursor.phasors.x_x.phase,2);
%                 phaseY.(groups{i})(:,k,j) = mean(data.(groups{i}){j}.(blocks{i,k}).cursor.phasors.y_y.phase,2);
%             end
%         end

        for j = 1:Nsubj
            for k = 1:3
                if i == 3 && j == 4 && k == 1
                    phaseX.(groups{i})(:,:,k,j) = [NaN(Nfreq,1) data.(groups{i}){j}.(blocks{i,k}).cursor.phasors.x_x.phase];
                    phaseY.(groups{i})(:,:,k,j) = [NaN(Nfreq,1) data.(groups{i}){j}.(blocks{i,k}).cursor.phasors.y_y.phase];
                else
                    phaseX.(groups{i})(:,:,k,j) = data.(groups{i}){j}.(blocks{i,k}).cursor.phasors.x_x.phase;
                    phaseY.(groups{i})(:,:,k,j) = data.(groups{i}){j}.(blocks{i,k}).cursor.phasors.y_y.phase;
                end
            end
        end
        
        delay.(groups{i}).x = -1000*phaseX.(groups{i})./(repmat(f_x,[1 5 3 Nsubj])*2*pi);
        delay.(groups{i}).y = -1000*phaseY.(groups{i})./(repmat(f_y,[1 5 3 Nsubj])*2*pi);
        
        delayMu.(groups{i}).x = mean(delay.(groups{i}).x,4,'omitnan');
        delayMu.(groups{i}).y = mean(delay.(groups{i}).y,4,'omitnan');
        
        delaySE.(groups{i}).x = std(delay.(groups{i}).x,[],4,'omitnan');
        delaySE.(groups{i}).y = std(delay.(groups{i}).y,[],4,'omitnan');
    end
    
    figure(1); clf
    for i = 1:length(groups)
        subplot(1,3,i); hold on
        for j = 1:Nblock
            for k = 1:Nfreq
                s = shadedErrorBar((j-1)*5 + (1:5),delayMu.(groups{i}).x(k,:,j),delaySE.(groups{i}).x(k,:,j));
                editErrorBar(s,col(k,:),1);
            end
        end
        
        title(names{i})
        xticks(1:5:11)
        xticklabels({'Baseline','Late','Flip'})
        ylim([-500 3000])
        if i == 1
            ylabel('Lag (ms)')
        elseif i == 2
            xlabel('Blocks')
        end
    end
end