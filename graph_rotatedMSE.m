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
            clear target cursor
            target(:,1) = data.(groups{p}).(subj{k}).(block_name{gblocks(j)}).target.x_pos;
            target(:,2) = data.(groups{p}).(subj{k}).(block_name{gblocks(j)}).target.y_pos;
            cursor(:,1) = data.(groups{p}).(subj{k}).(block_name{gblocks(j)}).cursor.x_pos;
            cursor(:,2) = data.(groups{p}).(subj{k}).(block_name{gblocks(j)}).cursor.y_pos;
            
            % bandpass signals if desired
%             target2 = bandpass(target2,bpFreqs,130.004)';
%             cursor = bandpass(cursor,bpFreqs,130.004)';
            
            % ensure that rows are x vs y and columns are points in trajectory
            if size(cursor,1) > size(cursor,2)
                cursor = cursor';
                target = target';
            end

            % rotate the cursor trajectories relative to the target
            for i = 1:length(ang)
                R = rotz(ang(i));
                R = R(1:2,1:2);
                Rcursor = R*cursor;
                
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
    minMSE_vec = reshape(minMSE,[numel(minMSE) 1]); % vectorize minMSE
    minAng = ang(idx); % rotation angle at which MSE is minimized
    minAng_avg = mean(minAng,2); % rotation angle averaged across participants
    minAng_vec = reshape(minAng,[numel(minAng) 1]); % vectorize minAng
    
    % calculate covariance matrix of minimum MSE vs rotation angle for confidence ellipse
    for i = 1:length(gblocks)
        covar(:,:,i) = cov(minAng(i,:),minMSE(i,:));
    end

    % plot distributions of minimum MSE vs rotation angle
    figure; hold on
%     colors = [1 5 2 3];
    for i = 1:length(gblocks)
        a = error_ellipse(covar(:,:,i),[minAng_avg(i),minMSE_avg(i)],'conf',0.5); % confidence ellipse
        a.Color = col(i,:);
    end
    scatter(minAng_avg,minMSE_avg,60,col(1:length(gblocks),:),'filled') % average minimum MSE
    scatter(minAng_vec,minMSE_vec,25,repmat(col(1:length(gblocks),:),[Nsubj 1]),'filled','MarkerFaceAlpha',0.5); % minimum MSE for each participant
%     scatter(minAng_avg,minMSE_avg,60,col(colors,:),'filled') % average minimum MSE
%     scatter(minAng_vec,minMSE_vec,25,repmat(col(colors,:),[Nsubj 1]),'filled','MarkerFaceAlpha',0.5); % minimum MSE for each participant
%     plot([-90 180],[nothing nothing],'--k','LineWidth',1) % MSE of doing nothing
    xticks(-60:30:180)
    yticks(0.002:0.0005:0.004)
    xlabel(['Rotation Angle (',char(176),')'])
    ylabel('Mean-Squared Error (m^2)')
%     axis([-45 110 .0017 .0035])
    set(gca,'TickDir','out')
    box off
    legend(graph_name{gblocks},'Location','Northwest')
    set(gcf,'Renderer','painters')

    se = std(minAng')/sqrt(Nsubj);
    minAng_after = minAng([1 end],:); % rotation angle from baseline and post-learning
end
end