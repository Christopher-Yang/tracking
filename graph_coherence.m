function graph_coherence(data)
% plots the coherence between target and hand movement

% set variables for plotting
output = 'Rhand';
groups = {'rot','mir'};
block_name = {'baseline','pert1','pert2','pert3','pert4','post'};
Ngroup = length(groups);
Nsubj = length(data.(groups{1}));
Ntrials = length(data.(groups{1}){1}.(block_name{1}).MSE);
Nblock = length(block_name);
Nfreq = length(data.(groups{1}){1}.(block_name{1}).freqX);
f_x = data.(groups{1}){1}.(block_name{1}).freqX;
f_y = data.(groups{1}){1}.(block_name{1}).freqY;
sorted_freqs = sort([f_x f_y]);
SR.x_all = NaN(Ntrials,Nfreq,Nblock,Nsubj,Ngroup);
SR.y_all = NaN(Ntrials,Nfreq,Nblock,Nsubj,Ngroup);

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
                coh = mscohere([outX(:,k) outY(:,k)],[targetX(:,k) ...
                    targetY(:,k)],blackmanharris(round(N/5)),[] ...
                    ,sorted_freqs,130.004,'mimo')';
                cohx(k,:) = coh(1,:); % store target-hand x coherence
                cohy(k,:) = coh(2,:); % store target-hand y coherence
            end
            
            % average across trials
            SR.x_all(:,:,j,i,p) = cohx(:,1:2:end); 
            SR.y_all(:,:,j,i,p) = cohy(:,2:2:end);
        end
    end
end

% mean and standard error of stimulus-response coherence across subjects
SR.x = squeeze(mean(SR.x_all,4));
SR.y = squeeze(mean(SR.y_all,4));
SR.xSE = squeeze(std(SR.x_all,[],4)/sqrt(Nsubj)); 
SR.ySE = squeeze(std(SR.y_all,[],4)/sqrt(Nsubj));

% generate Figure 4B and Figure 4-supplement 1B
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
end
