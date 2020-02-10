function graph_coherence(data)
% plots the coherence between target and hand movement

% set variables for plotting
output = 'Rhand';
gblocks = [1 2 5 6];
groups = {'rot','mir'};
block_name = {'baseline','pert1','pert2','pert3','pert4','post'};
graph_name = {'Baseline','Early','Train2','Train3','Late','Post'};
Nsubj = length(data.(groups{1}));
Ntrials = length(data.(groups{1}){1}.(block_name{1}).MSE);
Nblock = length(block_name);
f_x = data.(groups{1}){1}.(block_name{1}).freqX;
f_y = data.(groups{1}){1}.(block_name{1}).freqY;
sorted_freqs = sort([f_x f_y]);

% set line colors
col = lines;
col = col(1:7,:);

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
            
            % calculate coherence
            for k = 1:Ntrials
                coh = mscohere([outX(:,k) outY(:,k)],[targetX(:,k) ...
                    targetY(:,k)],blackmanharris(round(N/5)),[] ...
                    ,sorted_freqs,130.004,'mimo')';
                cohxx(k,:) = coh(1,:); % store target-hand x coherence
                cohyy(k,:) = coh(2,:); % store target-hand y coherence
            end
            
            % average across trials
            SR.x_all(j,:,i,p) = mean(cohxx(:,1:2:end)); 
            SR.y_all(j,:,i,p) = mean(cohyy(:,2:2:end));
        end
    end
end

% average coherence across subjects
SR.x = squeeze(mean(SR.x_all,3));
SR.y = squeeze(mean(SR.y_all,3));

% standard error of coherence across subjects
SR.xSE = squeeze(std(SR.x_all,[],3)/sqrt(Nsubj)); 
SR.ySE = squeeze(std(SR.y_all,[],3)/sqrt(Nsubj));

% generate Figure 4B and S2B
figure(15); clf
for j = 1:2
    subplot(2,2,j); hold on
    for i = 1:length(gblocks)
        s = shadedErrorBar(f_x,SR.x(gblocks(i),:,j)...
            ,SR.xSE(gblocks(i),:,j),'lineProps','-o');
        editErrorBar(s,col(i,:),1);
    end
    if j == 1
        title('Rotation')
    else
        title('Mirror-Reversal')
    end
    ylabel([output,' SR_X coherence'])
    yticks(0:0.2:1)
    axis([0 2.3 0 1])
    
    subplot(2,2,j+2); hold on
    for i = 1:length(gblocks)
        s = shadedErrorBar(f_y,SR.y(gblocks(i),:,j)...
            ,SR.ySE(gblocks(i),:,j),'lineProps','-o');
        editErrorBar(s,col(i,:),1);
    end
    xlabel('Frequency (Hz)')
    ylabel([output,' SR_Y coherence'])
    yticks(0:0.2:1)
    axis([0 2.3 0 1])
end
legend(graph_name(gblocks),'Location','southeast')
end
