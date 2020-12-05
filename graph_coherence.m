% function graph_coherence(data)
% plots the coherence between target and hand movement
data = data1;
% set variables for plotting
output = 'Rhand';
gblocks = [1 2 5 6];
groups = {'rot','mir'};
block_name = {'baseline','pert1','pert2','pert3','pert4','post'};
graph_name = {'Baseline','Early','Train2','Train3','Late','Post'};
Ngroup = length(groups);
Nsubj = length(data.(groups{1}));
Ntrials = length(data.(groups{1}){1}.(block_name{1}).MSE);
Nblock = length(block_name);
Nfreq = length(data.(groups{1}){1}.(block_name{1}).freqX);
f_x = data.(groups{1}){1}.(block_name{1}).freqX;
f_y = data.(groups{1}){1}.(block_name{1}).freqY;
sorted_freqs = sort([f_x f_y]);
% SR.x_all = NaN(Ntrials,Nfreq,Nblock,Nsubj,Ngroup);
% SR.y_all = NaN(Ntrials,Nfreq,Nblock,Nsubj,Ngroup);
RR.x_all = NaN(Nfreq,Nblock,Nsubj,Ngroup);
RR.y_all = NaN(Nfreq,Nblock,Nsubj,Ngroup);


% set line colors
col = copper;
col = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

% compute target-hand coherence (SR: stimulus-response)
for p = 1:length(groups)
    for i = 1:Nsubj
        for j = 1:Nblock
            
            % target position data in all trials
            targetX = data.(groups{p}){i}.(block_name{j}).target.x_pos; 
            targetY = data.(groups{p}){i}.(block_name{j}).target.y_pos;
            
            % hand position data in all trials
            outX = data.(groups{p}){i}.(block_name{j}).(output).x_pos; 
            outY = data.(groups{p}){i}.(block_name{j}).(output).y_pos;
            N = size(targetX,1);
            
            cohx = NaN(Ntrials, length(sorted_freqs));
            cohy = NaN(Ntrials, length(sorted_freqs));
            
            % calculate coherence
            for k = 1:Ntrials
%                 coh = mscohere([outX(:,k) outY(:,k)],[targetX(:,k) ...
%                     targetY(:,k)],blackmanharris(round(N/5)),[] ...
%                     ,sorted_freqs,130.004,'mimo')';
                coh = mscohere([outX(:,k) outY(:,k)],[targetX(:,k) ...
                    targetY(:,k)],rectwin(round(N/5)),[] ...
                    ,sorted_freqs,130.004,'mimo')';
                cohx(k,:) = coh(1,:); % store target-hand x coherence
                cohy(k,:) = coh(2,:); % store target-hand y coherence
            end
            
            % average across trials
            SR.x_all(:,:,j,i,p) = cohx(:,1:2:end); 
            SR.y_all(:,:,j,i,p) = cohy(:,2:2:end);
            
            FLAG = 1;
            k1 = 1;
            k2 = Ntrials;
            k3 = 1;
            cohx = NaN(nchoosek(4,2), length(sorted_freqs));
            cohy = NaN(nchoosek(4,2), length(sorted_freqs));
            while FLAG
%                 coh = mscohere([outX(:,k1) outY(:,k1)],[outX(:,k2) ...
%                     outY(:,k2)],blackmanharris(round(N/5)),[] ...
%                     ,sorted_freqs,130.004,'mimo')';
                coh = mscohere([outX(:,k1) outY(:,k1)],[outX(:,k2) ...
                    outY(:,k2)],rectwin(round(N/5)),[] ...
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
            RR.x_all(:,j,i,p) = mean(sqrt(cohx(:,1:2:end)),1);
            RR.y_all(:,j,i,p) = mean(sqrt(cohy(:,2:2:end)),1);
            
%             FLAG = 1;
%             k1 = Ntrials/2+1;
%             k2 = Ntrials;
%             k3 = 1;
%             cohx = NaN(nchoosek(4,2), length(sorted_freqs));
%             cohy = NaN(nchoosek(4,2), length(sorted_freqs));
%             while FLAG
%                 coh = mscohere([outX(:,k1) outY(:,k1)],[outX(:,k2) ...
%                     outY(:,k2)],blackmanharris(round(N/5)),[] ...
%                     ,sorted_freqs,130.004,'mimo')';
%                 cohx(k3,:) = coh(1,:);
%                 cohy(k3,:) = coh(2,:);
%                 k1 = k1 + 1;
%                 if k1 == k2
%                     k1 = Ntrials/2+1;
%                     k2 = k2 - 1;
%                 end
%                 if k2 == Ntrials/2+1
%                     FLAG = 0;
%                 end
%                 k3 = k3 + 1;
%             end
%             RR.x_all(2,:,j,i,p) = mean(sqrt(cohx(:,1:2:end)),1);
%             RR.y_all(2,:,j,i,p) = mean(sqrt(cohy(:,2:2:end)),1);
        end
    end
end

% mean and standard error of stimulus-response coherence across subjects
SR.x = squeeze(mean(SR.x_all,4));
SR.y = squeeze(mean(SR.y_all,4));
SR.xSE = squeeze(std(SR.x_all,[],4)/sqrt(Nsubj)); 
SR.ySE = squeeze(std(SR.y_all,[],4)/sqrt(Nsubj));

% mean and standard error of response-response coherence across subjects
% RR.x = squeeze(mean(RR.x_all,4));
% RR.y = squeeze(mean(RR.y_all,4));
% RR.xSE = squeeze(std(RR.x_all,[],4)/sqrt(Nsubj));
% RR.ySE = squeeze(std(RR.y_all,[],4)/sqrt(Nsubj));

% compute response that can be explained by nonlinear but not linear model
% SR.xBin = [mean(SR.x_all(1:4,:,:,:,:),1); mean(SR.x_all(5:8,:,:,:,:),1)];
% SR.yBin = [mean(SR.y_all(1:4,:,:,:,:),1); mean(SR.y_all(5:8,:,:,:,:),1)];
SR.xBin = squeeze(mean(SR.x_all,1));
SR.yBin = squeeze(mean(SR.y_all,1));
NL.x_all = RR.x_all - SR.xBin;
NL.y_all = RR.y_all - SR.yBin;
NL.x_all = NL.x_all(:,[1 2 5 6],:,:);
NL.y_all = NL.y_all(:,[1 2 5 6],:,:);
NL.x = squeeze(mean(NL.x_all,3));
NL.y = squeeze(mean(NL.y_all,3));
NL.xSE = squeeze(std(NL.x_all,[],3)/sqrt(Nsubj));
NL.ySE = squeeze(std(NL.y_all,[],3)/sqrt(Nsubj));

%%

% generate Figure 4B and S2B
figure(15); clf
for j = 1:2
    subplot(2,2,j); hold on
    rectangle('Position',[Ntrials+1 -0.5 4*Ntrials-1 2],'FaceColor',...
            [0 0 0 0.1],'EdgeColor','none');
    for i = 1:Nblock
        for k = 1:Nfreq
            plotIdx = Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials;
            s = shadedErrorBar(plotIdx,SR.x(:,k,i,j),SR.xSE(:,k,i,j));
            editErrorBar(s,col(k,:),1.5);
        end
    end
    if j == 1
        title('Rotation')
    else
        title('Mirror-Reversal')
    end
    set(gca,'TickDir','out')
    ylabel([output,' SR_X coherence'])
    xticks(1:8:41)
    yticks(0:0.25:1)
    axis([1 Ntrials*Nblock 0 1])
    
    subplot(2,2,j+2); hold on
    rectangle('Position',[Ntrials+1 -0.5 4*Ntrials-1 2],'FaceColor',...
            [0 0 0 0.1],'EdgeColor','none');
    for i = 1:Nblock
        for k = 1:Nfreq
            plotIdx = Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials;
            s = shadedErrorBar(plotIdx,SR.y(:,k,i,j),SR.ySE(:,k,i,j));
            editErrorBar(s,col(k,:),1.5);
        end
    end
    set(gca,'TickDir','out')
    xlabel('Trial Number')
    ylabel([output,' SR_Y coherence'])
    xticks(1:8:41) 
    yticks(0:0.25:1)
    axis([1 Ntrials*Nblock 0 1])
end

%%
interval = 8;
figure(16); clf
for j = 1:2
    subplot(2,2,j); hold on
%     rectangle('Position',[Ntrials+1 -0.5 4*Ntrials-1 2],'FaceColor',...
%             [0 0 0 0.1],'EdgeColor','none');
    plot([0 100],[0 0],'k')
%     for i = 1:Nblock
        for k = 1:Nfreq
%             plotIdx = 8*(i-1)*[1 1]+[1 5];
%             s = shadedErrorBar(plotIdx,NL.x(:,k,i,j),NL.xSE(:,k,i,j));
%             editErrorBar(s,col(k,:),1.5);
            plot((0:interval:3*interval)+k,NL.x(k,:,j),'.','Color',col(k,:),'MarkerSize',25)
            plot((0:interval:3*interval)+k,squeeze(NL.x_all(k,:,:,j)),'.','Color',col(k,:),'MarkerSize',10)
        end
%     end
    if j == 1
        title('Rotation')
    else
        title('Mirror-Reversal')
    end
    set(gca,'TickDir','out')
    ylabel('Nonlinear - linear coherence (x)')
    xticks(4:interval:64)
    xticklabels({'Baseline','Early','Late','Post'})
    yticks(-1:0.25:1)
    axis([0 32 -0.25 0.5])
    
    subplot(2,2,j+2); hold on
%     rectangle('Position',[Ntrials+1 -0.5 4*Ntrials-1 2],'FaceColor',...
%             [0 0 0 0.1],'EdgeColor','none');
    plot([0 100],[0 0],'k')
%     for i = 1:Nblock
        for k = 1:Nfreq
%             plotIdx = 8*(i-1)*[1 1]+[1 5];
%             s = shadedErrorBar(plotIdx,NL.y(:,k,i,j),NL.ySE(:,k,i,j));
%             editErrorBar(s,col(k,:),1.5);
            plot((0:interval:3*interval)+k,NL.y(k,:,j),'.','Color',col(k,:),'MarkerSize',25)
            plot((0:interval:3*interval)+k,squeeze(NL.y_all(k,:,:,j)),'.','Color',col(k,:),'MarkerSize',10)
        end
%     end
    set(gca,'TickDir','out')
    ylabel('Nonlinear - linear coherence (y)')
    xticks(4:interval:64)
    xticklabels({'Baseline','Early','Late','Post'})
    yticks(-0.25:0.25:0.5)
    axis([0 32 -0.25 0.5])
end
% end
