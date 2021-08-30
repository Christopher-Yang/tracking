function graph_coherence(data, block_name, graph_name)
% plots the coherence between target and hand movement

% set variables for plotting
groups = {'day2', 'day5', 'day10'}; % group names
group_names = {'2-day', '5-day', '10-day'};
Ngroup = length(groups); % number of groups
allSubj = [13 14 5]; % number of subjects

a = data.(groups{1}){1}.B1_baseline;
Ntrials = length(a.MSE); % number of trials
Nfreq = length(a.freqX); % number of frequencies
f_x = a.freqX; % x frequencies
f_y = a.freqY; % y frequencies
sorted_freqs = sort([f_x f_y]); % frequencies sorted in ascending order

% set line colors
col = copper;
col = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

% compute target-hand coherence (SR: stimulus-response coherence)
for p = 1:Ngroup
    
    Nblock = length(block_name.(groups{p})); % number of blocks

    dark{p} = find(contains(graph_name.(groups{p}),'D)'));
    flip{p} = find(contains(graph_name.(groups{p}),'(F)'));
    special{p} = find(contains(graph_name.(groups{p}),'('));
    normal{p} = 1:Nblock;
    normal{p}(special{p}) = [];
    
    % preallocate variables to store coherence
    SR{p}.x_all = NaN(Ntrials,Nfreq,Nblock,allSubj(p));
    SR{p}.y_all = NaN(Ntrials,Nfreq,Nblock,allSubj(p));
    RR{p}.x_all = NaN(Nfreq,Nblock,allSubj(p));
    RR{p}.y_all = NaN(Nfreq,Nblock,allSubj(p));

    for i = 1:allSubj(p)
        for j = 1:Nblock
            
            a = data.(groups{p}){i}.(block_name.(groups{p}){j});
            
            % target position data in all trials
            targetX = a.target.x_pos; 
            targetY = a.target.y_pos;
            
            % hand position data in all trials
            outX = a.cursor.x_pos; 
            outY = a.cursor.y_pos;
            N = size(targetX,1);
            
            % handles missing trial
            if p == 3 && i == 4 && j == 1
                Ntrials = 4;
            elseif p == 1 && i == 1 && j == 2
                Ntrials = 4;
            else
                Ntrials = 5;
            end
            
            % preallocate variables
            cohx = NaN(Ntrials, length(sorted_freqs));
            cohy = NaN(Ntrials, length(sorted_freqs));
            
            % calculate coherence
            for k = 1:Ntrials
                coh = mscohere([outX(:,k) outY(:,k)],[targetX(:,k) ...
                    targetY(:,k)],blackmanharris(round(N/5)),[] ...
                    ,sorted_freqs,130.004,'mimo')';
                cohx(k,:) = coh(1,:); % store target-hand x coherence
                cohy(k,:) = coh(2,:); % store target-hand y coherence
            end
            
            % average across trials
            SR{p}.x_all(1:Ntrials,:,j,i) = cohx(:,1:2:end); 
            SR{p}.y_all(1:Ntrials,:,j,i) = cohy(:,2:2:end);
            
            FLAG = 1;
            k1 = 1;
            k2 = Ntrials;
            k3 = 1;
            cohx = NaN(nchoosek(Ntrials,2), length(sorted_freqs));
            cohy = NaN(nchoosek(Ntrials,2), length(sorted_freqs));
            while FLAG
                coh = mscohere([outX(:,k1) outY(:,k1)],[outX(:,k2) ...
                    outY(:,k2)],blackmanharris(round(N/5)),[] ...
                    ,sorted_freqs,130.004,'mimo')';
                cohx(k3,:) = coh(1,:);
                cohy(k3,:) = coh(2,:);
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
            RR{p}.x_all(:,j,i) = mean(sqrt(cohx(:,1:2:end)),1);
            RR{p}.y_all(:,j,i) = mean(sqrt(cohy(:,2:2:end)),1);
        end
    end
    
    % mean and standard error of stimulus-response coherence across subjects
    SR{p}.x = squeeze(mean(SR{p}.x_all,4,'omitnan'));
    SR{p}.y = squeeze(mean(SR{p}.y_all,4,'omitnan'));
    SR{p}.xSE = squeeze(std(SR{p}.x_all,[],4,'omitnan')/sqrt(allSubj(p)));
    SR{p}.ySE = squeeze(std(SR{p}.y_all,[],4,'omitnan')/sqrt(allSubj(p)));
    
    RR{p}.x = squeeze(mean(RR{p}.x_all,3,'omitnan'));
    RR{p}.y = squeeze(mean(RR{p}.y_all,3,'omitnan'));
    RR{p}.xSE = squeeze(std(RR{p}.x_all,[],3,'omitnan')/sqrt(allSubj(p)));
    RR{p}.ySE = squeeze(std(RR{p}.y_all,[],3,'omitnan')/sqrt(allSubj(p)));
    
    x = squeeze(mean(SR{p}.x_all,1,'omitnan'));
    y = squeeze(mean(SR{p}.y_all,1,'omitnan'));
    
    NL{p}.x_all = RR{p}.x_all - x;
    NL{p}.y_all = RR{p}.y_all - y;
    
    NL{p}.x = squeeze(mean(NL{p}.x_all,3,'omitnan'));
    NL{p}.y = squeeze(mean(NL{p}.y_all,3,'omitnan'));
    NL{p}.xSE = squeeze(std(NL{p}.x_all,[],3,'omitnan')/sqrt(allSubj(p)));
    NL{p}.ySE = squeeze(std(NL{p}.y_all,[],3,'omitnan')/sqrt(allSubj(p)));
end

% generate figures
labels = {'Baseline', 'Early', 'Late', "Flip 1", "Flip 2"};
figure(13); clf
for j = 1:Ngroup
    Nblock = length(normal{j});
    
    % plot Figure 4B
    subplot(2,3,j); hold on
    for i = 1:Nblock+2
        if i <= Nblock
            for k = 1:Nfreq
                plotIdx = Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials;
                s = shadedErrorBar(plotIdx,SR{j}.x(:,k,normal{j}(i)),SR{j}.xSE(:,k,normal{j}(i)));
                editErrorBar(s,col(k,:),1.5); hold on
            end
        else
            for k = 1:Nfreq
                idx = i - Nblock;
                plotIdx = Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials;
                s = shadedErrorBar(plotIdx,SR{j}.x(:,k,flip{j}(idx)),SR{j}.xSE(:,k,flip{j}(idx)));
                editErrorBar(s,col(k,:),1.5); hold on
            end
        end
    end
    title(group_names{j})
    set(gca,'TickDir','out')
    if j == 1
        ylabel('Coherence for x-target freqs')
    end
    ticks = 1:5:(2+Nblock)*Ntrials;
    xticks(ticks([1 2 end-2 end-1 end]))
    xticklabels(labels)
    yticks(0:0.25:1)
    axis([1 Ntrials*(2+Nblock) 0 1])
    
    % plot Figure 4-supplement 1B
    subplot(2,3,j+3); hold on
    for i = 1:Nblock+2
        if i <= Nblock
            for k = 1:Nfreq
                plotIdx = Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials;
                s = shadedErrorBar(plotIdx,SR{j}.y(:,k,normal{j}(i)),SR{j}.ySE(:,k,normal{j}(i)));
                editErrorBar(s,col(k,:),1.5);
            end
        else
            for k = 1:Nfreq
                idx = i - Nblock;
                plotIdx = Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials;
                s = shadedErrorBar(plotIdx,SR{j}.y(:,k,flip{j}(idx)),SR{j}.ySE(:,k,flip{j}(idx)));
                editErrorBar(s,col(k,:),1.5); hold on
            end
        end
    end
    set(gca,'TickDir','out')
    if j == 1
        ylabel('Coherence for y-target freqs')
    end
    xticks(ticks([1 2 end-2 end-1 end]))
    xticklabels(labels)
    yticks(0:0.25:1)
    axis([1 Ntrials*(2+Nblock) 0 1])
end

labels = {'Baseline', 'Early', 'Late', 'Flip'};
figure(14); clf
for j = 1:Ngroup
    Nblock = length(dark{j});
    
    % plot Figure 4B
    subplot(2,3,j); hold on
    for i = 1:Nblock
        for k = 1:Nfreq
            plotIdx = Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials;
            s = shadedErrorBar(plotIdx,SR{j}.x(:,k,dark{j}(i)),SR{j}.xSE(:,k,dark{j}(i)));
            editErrorBar(s,col(k,:),1.5); hold on
        end
    end
    title(group_names{j})
    set(gca,'TickDir','out')
    if j == 1
        ylabel('Coherence for x-target freqs')
    end
    ticks = 1:5:Nblock*Ntrials;
    xticks(ticks([1 2 end-1 end]))
    xticklabels(labels)
    yticks(0:0.25:1)
    axis([1 Ntrials*Nblock 0 1])
    
    % plot Figure 4-supplement 1B
    subplot(2,3,j+3); hold on
    for i = 1:Nblock
        for k = 1:Nfreq
            plotIdx = Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials;
            s = shadedErrorBar(plotIdx,SR{j}.y(:,k,normal{j}(i)),SR{j}.ySE(:,k,normal{j}(i)));
            editErrorBar(s,col(k,:),1.5);
        end
    end
    set(gca,'TickDir','out')
    if j == 1
        ylabel('Coherence for y-target freqs')
    end
    xticks(ticks([1 2 end-1 end]))
    xticklabels(labels)
    yticks(0:0.25:1)
    axis([1 Ntrials*Nblock 0 1])
end

% figure(15); clf
% for j = 1:Ngroup
%     Nblock = length(block_name.groups({j}));
%     idx = 1:Nblock;
%     idx = idx([1 2 end-4]);
%     
%     % plot Figure 4B
%     subplot(2,3,j); hold on
%     for i = 1:3
%         for k = 1:Nfreq
%             plotIdx = Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials;
%             s = shadedErrorBar(plotIdx,SR{j}.x(:,k,idx(i)),SR{j}.xSE(:,k,idx(i)));
%             editErrorBar(s,col(k,:),1.5);
%         end
%     end
%     title(group_names{j})
%     set(gca,'TickDir','out')
%     if j == 1
%         ylabel('Coherence for x-target freqs')
%     end
%     xticks(1:5:11)
%     xticklabels({'Baseline','Early','Late'})
%     yticks(0:0.25:1)
%     axis([1 Ntrials*3 0 1])
%     
%     % plot Figure 4-supplement 1B
%     subplot(2,3,j+3); hold on
%     for i = 1:3
%         for k = 1:Nfreq
%             plotIdx = Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials;
%             s = shadedErrorBar(plotIdx,SR{j}.y(:,k,idx(i)),SR{j}.ySE(:,k,idx(i)));
%             editErrorBar(s,col(k,:),1.5);
%         end
%     end
%     set(gca,'TickDir','out')
%     if j == 1
%         ylabel('Coherence for y-target freqs')
%     end
%     xticks(1:5:11)
%     xticklabels({'Baseline','Early','Late'})
%     yticks(0:0.25:1)
%     axis([1 Ntrials*3 0 1])
% end
% 
% figure(16); clf
% for j = 1:Ngroup
%     % plot Figure 4B
%     subplot(2,3,j); hold on
%     for k = 1:Nfreq
%         errorbar(1:Nblock,NL{j}.x(k,:),NL{j}.xSE(k,:),'.','Color',col(k,:),'MarkerSize',20);
%     end
%     title(group_names{j})
%     set(gca,'TickDir','out')
%     if j == 1
%         ylabel('Coherence for x-target freqs')
%     end
%     xticks(1:Nblock)
%     xticklabels(labels)
%     yticks(0:0.25:1)
%     axis([0.5 Nblock+0.5 0 1])
%     
%     % plot Figure 4-supplement 1B
%     subplot(2,3,j+3); hold on
%     for k = 1:Nfreq
%         errorbar(1:Nblock,NL{j}.y(k,:),NL{j}.ySE(k,:),'.','Color',col(k,:),'MarkerSize',20);
%     end
%     set(gca,'TickDir','out')
%     if j == 1
%         ylabel('Coherence for y-target freqs')
%     end
%     xticks(1:Nblock)
%     xticklabels(labels)
%     yticks(0:0.25:1)
%     axis([0.5 Nblock+0.5 0 1])
% end

end