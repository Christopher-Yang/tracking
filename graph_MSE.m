% graph_MSE plots the mean-squared error between the target and cursor
% positions for each trial. 
% 
%   data: structure containing all data
%   block_name: name of blocks
%   blockType: identifies whether the block has visual feedback (=1), no 
%       visual feedback (=2), or the flipped mapping (=3)

function graph_MSE(data, block_name, blockType)

% set variables for plotting
groups = {'day2','day5','day10'};
Ntrials = length(data.(groups{1}){1}.(block_name.(groups{1}){1}).MSE);

% set color maps
col = [180 180 0
       0 191 255
       255 99 71]./255;

% store data in new variables
for i = 1:length(groups) % iterate over groups of subjects
    Nsubj = length(data.(groups{i}));
    normal{i} = find(blockType.(groups{i}) == 1); % with visual feedback
    dark{i} = find(blockType.(groups{i}) == 2); % without visual feedback
    
    Nblocks = length(fieldnames(data.(groups{i}){1}));
    
    % preallocate variables
    MSE_full = NaN(Ntrials, Nblocks, Nsubj);

    
    % store data in variables
    for j = 1:Nsubj % iterate over subjects
        x = NaN(Ntrials, Nblocks);
        for k = 1:Nblocks % iterate over blocks
            
            if (i == 1 && j == 1 && k == 2) || (i == 3 && j == 4 && k == 1)
                x(:,k) = [NaN data.(groups{i}){j}.(block_name.(groups{i}){k}).MSE];
            else
                x(:,k) = data.(groups{i}){j}.(block_name.(groups{i}){k}).MSE; % store MSE from each trial in x
            end
        end
        MSE_full(:,:,j) = x; % put the individual trial MSEs from x into MSE_full
    end
    MSE.(groups{i}) = permute(MSE_full, [1 3 2]);
end

offset = [0 20 55];

%% Figure 3B
f = figure(2); clf; hold on
set(f,'Position',[200 200 380 100]);
for i = 1:length(groups)
    Nblocks = length(normal{i});
    
    plot([offset(i) offset(i)]+1,[0 110],'k','LineWidth',0.5)
    plot([offset(i)+1 offset(i)+Ntrials*Nblocks],[10 10],'k','LineWidth',0.5)
    
    for j = 1:Nblocks
        plotIdx = ((j-1)*Ntrials+1:(j-1)*Ntrials+Ntrials) + offset(i);
        
        if j > 2
            plot([plotIdx(1)-0.5 plotIdx(1)-0.5],[0 110],'Color',[0.8 0.8 0.8])
        end
        
        plot(plotIdx,MSE.(groups{i})(:,:,normal{i}(j)),'Color',[col(i,:) 0.5],'LineWidth',0.3)
        plot(plotIdx,mean(MSE.(groups{i})(:,:,normal{i}(j)),2),'Color',col(i,:),'LineWidth',2)
    end
end
set(gca,'TickDir','out','Xcolor','none')
axis([1 Ntrials*Nblocks+offset(3) 10 110])
yticks(10:20:110)
ylabel('Mean squared error (cm)')

% save figure for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/MSE','-dpdf','-painters')

%% Supplementary Figure 2A
f = figure(3); clf; hold on
set(f,'Position',[200 200 380 100]);
for i = 1:length(groups)
    Nblocks = length(dark{i})-1;
    
    plot([offset(i) offset(i)]+1,[0 170],'k','LineWidth',0.5)
    plot([offset(i)+1 offset(i)+Ntrials*Nblocks],[20 20],'k','LineWidth',0.5)
    
    for j = 1:Nblocks
        plotIdx = ((j-1)*Ntrials+1:(j-1)*Ntrials+Ntrials) + offset(i);
        
        if j > 2
            plot([plotIdx(1)-0.5 plotIdx(1)-0.5],[0 170],'Color',[0.8 0.8 0.8])
        end
        
        plot(plotIdx,MSE.(groups{i})(:,:,dark{i}(j)),'Color',[col(i,:) 0.5],'LineWidth',0.3)
        plot(plotIdx,mean(MSE.(groups{i})(:,:,dark{i}(j)),2),'Color',col(i,:),'LineWidth',2)
    end
end
set(gca,'TickDir','out','Xcolor','none')
axis([1 Ntrials*Nblocks+offset(3) 20 170])
yticks(20:50:170)
ylabel('Mean squared error (cm)')

% save figure for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/MSE_dark','-dpdf','-painters')
