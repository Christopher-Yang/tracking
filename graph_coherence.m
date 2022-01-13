% graph_coherence plots the coherence between target and cursor movement
% 
%   data: structure containing all data
%   block_name: name of blocks
%   blockType: identifies whether the block has visual feedback (=1), no 
%       visual feedback (=2), or the flipped mapping (=3)

function graph_coherence(data, block_name, blockType)

% set variables for plotting
groups = {'day2', 'day5', 'day10'}; % group names
group_names = {'2-day', '5-day', '10-day'};
Ngroup = length(groups); % number of groups
allSubj = [13 14 5]; % number of subjects

a = data.(groups{1}){1}.B1_baseline;
Ntrials = size(a.cursor.x_pos,2); % number of trials
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
    
    % identify trial types
    normal{p} = find(blockType.(groups{p}) == 1); % with visual feedback
    dark{p} = find(blockType.(groups{p}) == 2); % without visual feedback
    
    % preallocate variables to store coherence
    SR{p}.x_all = NaN(Ntrials,Nfreq,Nblock,allSubj(p));
    SR{p}.y_all = NaN(Ntrials,Nfreq,Nblock,allSubj(p));

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
            
            N = size(targetX,1); % number of data samples in one trial
            
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
            
        end
    end
    
    % mean and standard error of stimulus-response coherence across subjects
    SR{p}.x = squeeze(mean(SR{p}.x_all,4,'omitnan'));
    SR{p}.y = squeeze(mean(SR{p}.y_all,4,'omitnan'));
    SR{p}.xSE = squeeze(std(SR{p}.x_all,[],4,'omitnan')/sqrt(allSubj(p)));
    SR{p}.ySE = squeeze(std(SR{p}.y_all,[],4,'omitnan')/sqrt(allSubj(p)));
end

disp('Done')

%% plot data from trials with visual feedback
offset = [0 20 55];

% generate figures
f = figure(1); clf
set(f,'Position',[200 200 380 250]);
for j = 1:Ngroup
    Nblock = length(normal{j});
    
    % plot Figure 3B; x-axis coherence
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
    
    % y-axis coherence (not shown in paper)
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

% save figure for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/coherence_normal','-dpdf','-painters')

blockNum = [3 5; 3 11; 3 21];

y = [];
for k = 1:Ngroup
    for i = 1:2        
        t = mean(SR{k}.x_all(:,:,blockNum(k,i),:),1,'omitnan');
        y = [y; t(:)];
    end
end

g(1:156,1) = "2-day";
g(157:324,1) = "5-day";
g(325:384,1) = "10-day";
b([1:78 157:240 325:354],1) = "early";
b([79:156 241:324 355:384],1) = "late";
frequency = repmat((1:Nfreq)',[64 1]);
s1 = repmat(1:13,[Nfreq 2]);
s2 = repmat(14:27,[Nfreq 2]);
s3 = repmat(28:32,[Nfreq 2]);
subject = [s1(:); s2(:); s3(:)];
T = table(g, b, frequency, subject, y, 'VariableNames', {'group','block','frequency','subject','coherence'});
writetable(T,'C:/Users/Chris/Documents/R/habit/data/coherence.csv')

%% plot data from trials without visual feedback
f = figure(2); clf
set(f,'Position',[200 200 380 250]);
for j = 1:Ngroup
    Nblock = length(dark{j});
    
    % plot Supplementary Figure 2; x-axis coherence
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
    
    % y-axis coherence (not shown in paper)
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
        ylabel('SR_y')
        yticks(0:0.25:1)
        axis([1 Ntrials*Nblock+offset(3) 0 1])
    end
end

% save for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/coherence_dark','-dpdf','-painters')

end