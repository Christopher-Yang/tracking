function graph_rotatedMSE(data,groups,block_name,gblocks,graph_name)

rng(1)
col = lines;
col = col(1:7,:);
delt = 0.01;
t = 0:delt:10-delt;
nstep = length(t);
ang = -90:1:180; % angles to rotate the cursor trajectory
Nsubj = length(fieldnames(data.(groups{1})))-1;
% bpFreqs = [1.5 2]; % frequencies to bandpass trajectories; optional
% lag = 0; % lag cursor trajectories; optional
% Nlag = round(lag/(1/130.004)); 

for p = 1:length(groups)
    for k = 1:Nsubj
        subj = fieldnames(data.(groups{p}));
        subj(end) = [];
        for j = 1:length(gblocks)
            % add all trajectories into target
            clear target hand
            target(:,1) = data.(groups{p}).(subj{k}).(block_name{gblocks(j)}).target.x_pos;
            target(:,2) = data.(groups{p}).(subj{k}).(block_name{gblocks(j)}).target.y_pos;
            hand(:,1) = data.(groups{p}).(subj{k}).(block_name{gblocks(j)}).Rhand.x_pos;
            hand(:,2) = data.(groups{p}).(subj{k}).(block_name{gblocks(j)}).Rhand.y_pos;
            
            % bandpass signals if desired
%             target2 = bandpass(target2,bpFreqs,130.004)';
%             cursor = bandpass(cursor,bpFreqs,130.004)';
            
            % ensure that rows are x vs y and columns are points in trajectory
            if size(hand,1) > size(hand,2)
                hand = hand';
                target = target';
            end

            % rotate the cursor trajectories relative to the target
            for i = 1:length(ang)
                R = rotz(ang(i));
                R = R(1:2,1:2);
                Rcursor = R*hand;
                
                % lag the target and cursor trajectories if desired; uncomment lag in line 10/11 
%                 if lag ~= 0
%                     target = target(:,1:end-Nlag);
%                     Rcursor = Rcursor(:,1:end-Nlag);
%                 end
                
                err = (Rcursor - target).^2;
                MSE(i,j,k) = mean(sum(err,1));
            end
        end
    end

    nothing = mean(sum(target.^2,1)); % MSE of doing nothing
    [minMSE,idx] = min(MSE); % minimum MSE for each block of tracking for each participant
    idx = squeeze(idx);
    minMSE = squeeze(minMSE);
    minMSE_avg = mean(minMSE,2); % minimium MSE averaged across participants
    minAng = ang(idx); % rotation angle at which MSE is minimized
    minAng_avg = mean(minAng,2); % rotation angle averaged across participants
    
    % calculate covariance matrix of minimum MSE vs rotation angle for confidence ellipse
    for i = 1:length(gblocks)
        covar(:,:,i) = cov(minAng(i,:),minMSE(i,:));
    end
    
    bl = [1 4];
    if p == 1
        minMSE_vec = reshape(minMSE,[numel(minMSE) 1]); % vectorize minMSE
        minAng_vec = reshape(minAng,[numel(minAng) 1]); % vectorize minAng
    else
        minMSE_vec = reshape(minMSE(bl,:),[numel(minMSE(bl,:)) 1]); % vectorize minMSE
        minAng_vec = reshape(minAng(bl,:),[numel(minAng(bl,:)) 1]); % vectorize minAng
        minMSE_avg = minMSE_avg(bl);
        minAng_avg = minAng_avg(bl);
        covar = covar(:,:,bl);
    end
    
    % plot distributions of minimum MSE vs rotation angle
    figure; hold on
    if p == 1
        for i = 1:length(gblocks)
            if i == 1
                plot(-1*ones(2,4),'.','MarkerSize',20)
            end
            hold on
            a = error_ellipse(covar(:,:,i),[minAng_avg(i),minMSE_avg(i)],'conf',0.5); % confidence ellipse
            a.Color = col(i,:);
            scatter(minAng_avg,minMSE_avg,60,col(1:length(gblocks),:),'filled') % average minimum MSE
            scatter(minAng_vec,minMSE_vec,25,repmat(col(1:length(gblocks),:),[Nsubj 1]),'filled','MarkerFaceAlpha',0.5); % minimum MSE for each participant
            title('Rotation')
            legend(graph_name{gblocks},'Location','Northwest')
        end
    else
        for i = 1:2
            if i == 1
                plot([-1 -1],'.','MarkerSize',20,'MarkerEdgeColor',col(1,:))
                hold on
                plot([-1 -1],'.','MarkerSize',20,'MarkerEdgeColor',col(4,:))
            end
            a = error_ellipse(covar(:,:,i),[minAng_avg(i),minMSE_avg(i)],'conf',0.5); % confidence ellipse
            a.Color = col(bl(i),:);
            scatter(minAng_avg,minMSE_avg,60,col(bl,:),'filled') % average minimum MSE
            scatter(minAng_vec,minMSE_vec,25,repmat(col(bl,:),[Nsubj 1]),'filled','MarkerFaceAlpha',0.5); % minimum MSE for each participant
            title('Mirror Reversal')
            legend(graph_name{gblocks(bl)},'Location','Northwest')
        end
    end
%     plot([-90 180],[nothing nothing],'--k','LineWidth',1) % MSE of doing nothing
    xticks(-90:30:90)
    yticks(0.002:0.0005:0.004)
    xlabel(['Rotation Angle (',char(176),')'])
    ylabel('Mean-Squared Error (m^2)')
    axis([-90 30 .0017 .0035])
    set(gca,'TickDir','out')
    box off
    set(gcf,'Renderer','painters')

    se = std(minAng')/sqrt(Nsubj);
    minAng_after = minAng([1 end],:); % rotation angle from baseline and post-learning
    
    % plot rotation angle distributions at different timepoints
    figure(2)
    if p == 1 % plot VMR data
        scatter(repmat((1:2)',[Nsubj 1]),reshape(minAng_after,[20 1]),25,repmat(col([1 4],:),[Nsubj 1]),'filled','MarkerFaceAlpha',0.5)
        hold on
        plot(repmat((1:2)',[1 10]),minAng_after,'Color',[0 0 0 0.5])
        errorbar(1,minAng_avg(1),2*se(1),'.','Color',col(1,:),'MarkerSize',30,'LineWidth',1)
        errorbar(2,minAng_avg(end),2*se(end),'.','Color',col(length(gblocks),:),'MarkerSize',30,'LineWidth',1)
        plot([0 9],[0 0],'--k','LineWidth',1)
    else % plot MR data
        scatter(repmat([3.5 4.5]',[Nsubj 1]),reshape(minAng_after,[numel(minAng_after) 1]),25,repmat(col([1 length(gblocks)],:),[Nsubj 1]),'filled','MarkerFaceAlpha',0.5)
        plot(repmat((3.5:4.5)',[1 Nsubj]),minAng_after,'Color',[0 0 0 0.5])
        errorbar(3.5,minAng_avg(1),2*se(1),'.','Color',col(1,:),'MarkerSize',30,'LineWidth',1)
        errorbar(4.5,minAng_avg(end),2*se(end),'.','Color',col(length(gblocks),:),'MarkerSize',30,'LineWidth',1)
    end
    xticks([1:2 3.5:4.5])
    axis([0.5 5 -50 30])
    xticklabels({'Baseline','Post','Baseline','Post'})
    yticks(-60:10:180)
    ylabel(['Rotation Angle (',char(176),')'])
    set(gca,'TickDir','out')
    box off
    
    minima(:,:,p) = minAng'; % save rotation angles in a new matrix for statistical analysis
end

% comp = [minima(:,[1 4],1); minima(:,[1 4],2)];
% comp2 = [repelem(0,10)'; repelem(1,10)'];
% within = {'time'};
% between = {'group'};
% [tbl,rm] = simple_mixed_anova(comp, comp2, within, between)

end