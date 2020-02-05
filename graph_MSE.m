function graph_MSE(data)

% set variables for plotting
groups = {'rot','mir'};
block_name = {'baseline','pert1','pert2','pert3','pert4','post'};
Nsubj = length(data.(groups{1}))-1;
Nblocks = length(block_name);
Ntrials = length(data.(groups{1}){1}.(block_name{1}).MSE);

% set color maps
col = lines;
col = [col(1:7,:) repelem(0.1,7)'];
col2 = [160 82 45
    46 139 87
    65 105 225]./255;

% store data in new variables
for i = 1:length(groups)
    % preallocate variables
    MSE_subj = NaN(Nblocks, Nsubj);
    MSE_full = NaN(Nblocks*Ntrials, Nsubj);
    
    % store data in variables
    for j = 1:Nsubj % iterate over subjects
        x = NaN(Nblocks*Ntrials, 1); 
        for k = 1:Nblocks % iterate over blocks
            MSE_subj(k,j) = mean(data.(groups{i}){j}.(block_name{k}).MSE); % average MSE within subject
            x((k-1)*Ntrials+1:(k-1)*Ntrials+Ntrials) = data.(groups{i}){j}.(block_name{k}).MSE; % store MSE from each trial in x
        end
        MSE_full(:,j) = x; % put the individual trial MSEs from x into MSE_full
    end
    MSE.(groups{i}).full = MSE_full;
    MSE.(groups{i}).full_se = std(MSE_full,[],2)/sqrt(size(MSE_full,2)); % compute standard error across subjects
end

% generate Figure 2C
figure(2); clf
for i = 1:6 % draw shaded error bar in 6 segments corresponding to 6 blocks
    s = shadedErrorBar(8*(i-1)+1:8*(i-1)+8,mean(MSE.rot.full(8*(i-1)+1:8*(i-1)+8,:),2),MSE.rot.full_se(8*(i-1)+1:8*(i-1)+8)); hold on;
    editErrorBar(s,col2(1,:),1);
    s = shadedErrorBar(8*(i-1)+1:8*(i-1)+8,mean(MSE.mir.full(8*(i-1)+1:8*(i-1)+8,:),2),MSE.mir.full_se(8*(i-1)+1:8*(i-1)+8));
    editErrorBar(s,col2(2,:),1);
end
rectangle('Position',[1 0 Ntrials 20],'FaceColor',col(1,:),'EdgeColor','none') % color different blocks with rectangles
rectangle('Position',[Ntrials+1 0 Ntrials 20],'FaceColor',col(2,:),'EdgeColor','none')
rectangle('Position',[Ntrials*4+1 0 Ntrials 20],'FaceColor',col(3,:),'EdgeColor','none')
rectangle('Position',[Ntrials*5+1 0 Ntrials 20],'FaceColor',col(4,:),'EdgeColor','none')
set(gca,'xtick',1:Ntrials:length(MSE.(groups{1}).full),'ytick',0:0.003:0.012,'TickDir','out')
axis([1 Ntrials*Nblocks 0 0.012])
xlabel('Tracking Cycle Number')
ylabel('Mean Squared Error (m^2)');
legend('Rotation','Mirror Reversal');
box off
legend boxoff
end