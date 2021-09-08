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

disp('Analyzing...')

% compute target-hand coherence (SR: stimulus-response coherence)
for p = 1:Ngroup
    
    disp(['   ' group_names{p} ' group'])
    
    Nblock = length(block_name.(groups{p})); % number of blocks

    dark{p} = find(contains(graph_name.(groups{p}),'(D)'));
    darkFlip{p} = find(contains(graph_name.(groups{p}),'(FD)'));
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
        
        disp(['      subject ' num2str(i)])
        
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

%%
offset = [0 20 55];

% generate figures
f = figure(1); clf
set(f,'Position',[200 200 380 250]);
for j = 1:Ngroup
    Nblock = length(normal{j});
    
    % plot Figure 4B
    subplot(2,1,1); hold on
    
    plot([offset(j) offset(j)]+1,[0 1],'k')
    plot([offset(j)+1 offset(j)+Ntrials*Nblock],[0 0],'k','LineWidth',0.5)
    for i = 1:Nblock
        plotIdx = (Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials) + offset(j);
        if i > 2
            plot([plotIdx(1)-0.5 plotIdx(1)-0.5],[0 1],'Color',[0.8 0.8 0.8])
        end
        for k = 1:Nfreq
            s = shadedErrorBar(plotIdx,SR{j}.x(:,k,normal{j}(i)),SR{j}.xSE(:,k,normal{j}(i)));
            editErrorBar(s,col(k,:),1.5); hold on
        end
    end
    if j == 3
        set(gca,'TickDir','out','Xcolor','none')
        ylabel('SR_x')
        yticks(0:0.25:1)
        axis([1 Ntrials*Nblock+offset(3) 0 1])
    end
    
    % plot Figure 4-supplement 1B
    subplot(2,1,2); hold on

    plot([offset(j) offset(j)]+1,[0 1],'k')
    plot([offset(j)+1 offset(j)+Ntrials*Nblock],[0 0],'k','LineWidth',0.5)
    for i = 1:Nblock
        plotIdx = (Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials) + offset(j);
        if i > 2
            plot([plotIdx(1)-0.5 plotIdx(1)-0.5],[0 1],'Color',[0.8 0.8 0.8])
        end
        for k = 1:Nfreq
            s = shadedErrorBar(plotIdx,SR{j}.y(:,k,normal{j}(i)),SR{j}.ySE(:,k,normal{j}(i)));
            editErrorBar(s,col(k,:),1.5);
        end
    end
    if j == 3
        set(gca,'TickDir','out','Xcolor','none')
        ylabel('SR_y')
        yticks(0:0.25:1)
        axis([1 Ntrials*Nblock+offset(3) 0 1])
    end
end

print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/coherence_normal','-dpdf','-painters')

%%
labels = {'Late',"Flip 1","Flip 2"};
f = figure(2); clf
set(f,'Position',[200 200 400 250]);
for j = 1:Ngroup
    Nblock = length(flip{j})+1;
    
    % plot Figure 4B
    subplot(2,3,j); hold on
    for i = 1:Nblock
        
        plotIdx = Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials;
        if i == 1
            block = normal{j}(end);
        else
            block = flip{j}(i-1);
        end
        
        for k = 1:Nfreq
            s = shadedErrorBar(plotIdx,SR{j}.x(:,k,block),SR{j}.xSE(:,k,block));
            editErrorBar(s,col(k,:),1.5); hold on
        end
    end
    title(group_names{j})
    set(gca,'TickDir','out')
    if j == 1
        ylabel('SR_x')
    end
    ticks = 1:5:Nblock*Ntrials;
    xticks(ticks([1 2 end]))
    xticklabels(labels)
    yticks(0:0.25:1)
    axis([1 Ntrials*Nblock 0 1])
    
    % plot Figure 4-supplement 1B
    subplot(2,3,j+3); hold on
    for i = 1:Nblock
        
        plotIdx = Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials;
        
        if i == 1
            block = normal{j}(end);
        else
            block = flip{j}(i-1);
        end
        
        for k = 1:Nfreq
            s = shadedErrorBar(plotIdx,SR{j}.y(:,k,block),SR{j}.ySE(:,k,block));
            editErrorBar(s,col(k,:),1.5);
        end
    end
    set(gca,'TickDir','out')
    if j == 1
        ylabel('SR_y')
    end
    xticks(ticks([1 2 end]))
    xticklabels(labels)
    yticks(0:0.25:1)
    axis([1 Ntrials*Nblock 0 1])
end

print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/coherence_flip','-dpdf','-painters')

%%
offset = [0 20 55];

f = figure(3); clf
set(f,'Position',[200 200 380 250]);
for j = 1:Ngroup
    Nblock = length(dark{j});
    
    % plot Figure 4B
    subplot(2,1,1); hold on

    plot([offset(j) offset(j)]+1,[0 1],'k')
    plot([offset(j)+1 offset(j)+Ntrials*Nblock],[0 0],'k','LineWidth',0.5)
    for i = 1:Nblock
        plotIdx = (Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials) + offset(j);
        if i > 2
            plot([plotIdx(1)-0.5 plotIdx(1)-0.5],[0 1],'Color',[0.8 0.8 0.8])
        end
        for k = 1:Nfreq
            s = shadedErrorBar(plotIdx,SR{j}.x(:,k,dark{j}(i)),SR{j}.xSE(:,k,dark{j}(i)));
            editErrorBar(s,col(k,:),1.5); hold on
        end
    end
    if j == 3
        set(gca,'TickDir','out','Xcolor','none')
        ylabel('SR_x')
        yticks(0:0.25:1)
        axis([1 Ntrials*Nblock+offset(3) 0 1])
    end
    
    % plot Figure 4-supplement 1B
    subplot(2,1,2); hold on
    
    plot([offset(j) offset(j)]+1,[0 1],'k')
    plot([offset(j)+1 offset(j)+Ntrials*Nblock],[0 0],'k','LineWidth',0.5)
    for i = 1:Nblock
        plotIdx = (Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials) + offset(j);
        if i > 2
            plot([plotIdx(1)-0.5 plotIdx(1)-0.5],[0 1],'Color',[0.8 0.8 0.8])
        end
        for k = 1:Nfreq
            s = shadedErrorBar(plotIdx,SR{j}.y(:,k,dark{j}(i)),SR{j}.ySE(:,k,dark{j}(i)));
            editErrorBar(s,col(k,:),1.5);
        end
    end
    if j == 3
        set(gca,'TickDir','out','Xcolor','none')
        ylabel('SR_x')
        yticks(0:0.25:1)
        axis([1 Ntrials*Nblock+offset(3) 0 1])
    end
end

print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/coherence_dark','-dpdf','-painters')

%%
labels = {'Late','Flip'};
f = figure(4); clf
set(f,'Position',[200 200 400 250]);
for j = 1:Ngroup
    Nblock = length(darkFlip{j})+1;
    
    % plot Figure 4B
    subplot(2,3,j); hold on
    for i = 1:Nblock
        
        plotIdx = Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials;
        if i == 1
            block = dark{j}(end);
        else
            block = darkFlip{j}(i-1);
        end
        
        for k = 1:Nfreq
            s = shadedErrorBar(plotIdx,SR{j}.x(:,k,block),SR{j}.xSE(:,k,block));
            editErrorBar(s,col(k,:),1.5); hold on
        end
    end
    title(group_names{j})
    set(gca,'TickDir','out')
    if j == 1
        ylabel('SR_x')
    end
    ticks = 1:5:Nblock*Ntrials;
    xticks(ticks)
    xticklabels(labels)
    yticks(0:0.25:1)
    axis([1 Ntrials*Nblock 0 1])
    
    % plot Figure 4-supplement 1B
    subplot(2,3,j+3); hold on
    for i = 1:Nblock
        
        plotIdx = Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials;
        if i == 1
            block = dark{j}(end);
        else
            block = darkFlip{j}(i-1);
        end
        
        for k = 1:Nfreq
            s = shadedErrorBar(plotIdx,SR{j}.y(:,k,block),SR{j}.ySE(:,k,block));
            editErrorBar(s,col(k,:),1.5);
        end
    end
    set(gca,'TickDir','out')
    if j == 1
        ylabel('SR_y')
    end
    xticks(ticks)
    xticklabels(labels)
    yticks(0:0.25:1)
    axis([1 Ntrials*Nblock 0 1])
end

print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/coherence_darkFlip','-dpdf','-painters')

%%

labels = {'Baseline', 'Early', 'Late'};
figure(5); clf
for j = 1:Ngroup
    Nblock = length(normal{j});
    
    % plot Figure 4B
    subplot(2,3,j); hold on
    for k = 1:Nfreq
        s = shadedErrorBar(1:Nblock,NL{j}.x(k,normal{j}),NL{j}.xSE(k,normal{j}));
        editErrorBar(s,col(k,:),1.5)
    end
    title(group_names{j})
    set(gca,'TickDir','out')
    if j == 1
        ylabel('NL_x')
    end
    ticks = 1:Nblock;
    xticks(ticks([1 2 end]))
    xticklabels(labels)
    yticks(0:0.25:1)
    axis([1 Nblock 0 1])
    
    % plot Figure 4-supplement 1B
    subplot(2,3,j+3); hold on
    for k = 1:Nfreq
        s = shadedErrorBar(1:Nblock,NL{j}.y(k,normal{j}),NL{j}.ySE(k,normal{j}));
        editErrorBar(s,col(k,:),1.5)
    end
    set(gca,'TickDir','out')
    if j == 1
        ylabel('NL_y')
    end
    xticks(ticks([1 2 end]))
    xticklabels(labels)
    yticks(0:0.25:1)
    axis([1 Nblock 0 1])
end

%%
labels = {'Late',"Flip 1","Flip 2"};

figure(6); clf
for j = 1:Ngroup
    Nblock = length(flip{j})+1;
    normal_end = normal{j}(end);
    
    % plot Figure 4B
    subplot(2,3,j); hold on
    for k = 1:Nfreq
        s = shadedErrorBar(1:Nblock,NL{j}.x(k,[normal_end flip{j}]),NL{j}.xSE(k,[normal_end flip{j}]));
        editErrorBar(s,col(k,:),1.5)
    end
    title(group_names{j})
    set(gca,'TickDir','out')
    if j == 1
        ylabel('NL_x')
    end
    xticks(1:Nblock)
    xticklabels(labels)
    yticks(0:0.25:1)
    axis([1 Nblock 0 1])
    
    % plot Figure 4-supplement 1B
    subplot(2,3,j+3); hold on
    for k = 1:Nfreq
        s = shadedErrorBar(1:Nblock,NL{j}.y(k,[normal_end flip{j}]),NL{j}.ySE(k,[normal_end flip{j}]));
        editErrorBar(s,col(k,:),1.5)
    end
    set(gca,'TickDir','out')
    if j == 1
        ylabel('NL_y')
    end
    xticks(1:Nblock)
    xticklabels(labels)
    yticks(0:0.25:1)
    axis([1 Nblock 0 1])
end

%%

labels = {'Baseline', 'Early', 'Late'};
figure(7); clf
for j = 1:Ngroup
    Nblock = length(dark{j});
    
    % plot Figure 4B
    subplot(2,3,j); hold on
    for k = 1:Nfreq
        s = shadedErrorBar(1:Nblock,NL{j}.x(k,dark{j}),NL{j}.xSE(k,dark{j}));
        editErrorBar(s,col(k,:),1.5)
    end
    title(group_names{j})
    set(gca,'TickDir','out')
    if j == 1
        ylabel('NL_x')
    end
    ticks = 1:Nblock;
    xticks(ticks([1 2 end]))
    xticklabels(labels)
    yticks(0:0.25:1)
    axis([1 Nblock 0 1])
    
    % plot Figure 4-supplement 1B
    subplot(2,3,j+3); hold on
    for k = 1:Nfreq
        s = shadedErrorBar(1:Nblock,NL{j}.x(k,dark{j}),NL{j}.xSE(k,dark{j}));
        editErrorBar(s,col(k,:),1.5)
    end
    set(gca,'TickDir','out')
    if j == 1
        ylabel('NL_y')
    end
    xticks(ticks([1 2 end]))
    xticklabels(labels)
    yticks(0:0.25:1)
    axis([1 Nblock 0 1])
end


end