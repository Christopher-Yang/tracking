function graph_MSE(data)
% plots the mean-squared error of every tracking trial

% set variables for plotting
groups = {'rot','mir'};
block_name = {'baseline','pert1','pert2','pert3','pert4','post'};
Nsubj = length(data.(groups{1}));
Nblocks = length(block_name);
Ntrials = length(data.(groups{1}){1}.(block_name{1}).MSE);

% set color maps
col = lines;
col = col(1:6,:);
col(5:6,:) = col(3:4,:);
col(3:4,:) = [0 0 0; 0 0 0];

% store data in new variables
for i = 1:length(groups) % iterate over groups of subjects
    % preallocate variables
    MSE_full = NaN(Ntrials, Nblocks, Nsubj);
    
    % store data in variables
    for j = 1:Nsubj % iterate over subjects
        x = NaN(Ntrials, Nblocks); 
        for k = 1:Nblocks % iterate over blocks
            x(:,k) = data.(groups{i}){j}.(block_name{k}).MSE; % store MSE from each trial in x
        end
        MSE_full(:,:,j) = x; % put the individual trial MSEs from x into MSE_full
    end
    MSE.(groups{i}) = permute(MSE_full, [1 3 2]);
end

% generate Figure 2C
figure(2); clf
subplot(1,2,1); hold on
rectangle('Position',[9 0 31 150],'FaceColor',[0 0 0 0.1],'EdgeColor','none');
for i = 1:6
    plot((i-1)*8+1:(i-1)*8+8,MSE.rot(:,:,i),'Color',[0 0 0 0.5],'LineWidth',0.2)
    plot((i-1)*8+1:(i-1)*8+8,mean(MSE.rot(:,:,i),2),'Color',col(i,:),'LineWidth',1.5)
end
title('Rotation')
xlabel('Trial Number (tracking)')
ylabel('Mean-squared error (cm^2)')
xticks(1:8:41)
axis([1 49 0 150])
set(gca,'TickDir','out')

subplot(1,2,2); hold on
rectangle('Position',[9 0 31 150],'FaceColor',[0 0 0 0.1],'EdgeColor','none');
for i = 1:6
    plot((i-1)*8+1:(i-1)*8+8,MSE.mir(:,:,i),'Color',[0 0 0 0.5],'LineWidth',0.2)
    plot((i-1)*8+1:(i-1)*8+8,mean(MSE.mir(:,:,i),2),'Color',col(i,:),'LineWidth',1.5)
end
title('Mirror Reversal')
xlabel('Trial Number (tracking)')
xticks(1:8:41)
axis([1 49 0 150])
set(gca,'TickDir','out')

