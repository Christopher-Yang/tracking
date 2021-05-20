function graph_coherence(data, output)
% plots the coherence between target and hand movement

% set variables for plotting
groups = {'day2', 'day5', 'day10'}; % group names
graph_names = {'2-day', '5-day', '10-day'};
block_name{1} = {'B1_baseline', 'B3', 'B9', 'B11_habit'}; % block names
block_name{2} = {'B1_baseline', 'B3', 'B24', 'B26_habit'}; % block names
block_name{3} = {'B1_baseline', 'B3', 'B49', 'B51_habit'}; % block names
Ngroup = length(groups); % number of groups
Nblock = length(block_name{1}); % number of blocks
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

% compute target-hand coherence (SR: stimulus-response coherence)
for p = 1:length(groups)
    
    % preallocate variables to store coherence
    SR{p}.x_all = NaN(Ntrials,Nfreq,Nblock,allSubj(p));
    SR{p}.y_all = NaN(Ntrials,Nfreq,Nblock,allSubj(p));
    RR{p}.x_all = NaN(Nfreq,Nblock,allSubj(p));
    RR{p}.y_all = NaN(Nfreq,Nblock,allSubj(p));

    for i = 1:allSubj(p)
        for j = 1:Nblock
            
            a = data.(groups{p}){i}.(block_name{p}{j});
            
            % target position data in all trials
            targetX = a.target.x_pos; 
            targetY = a.target.y_pos;
            
            % hand position data in all trials
            outX = a.(output).x_pos; 
            outY = a.(output).y_pos;
            N = size(targetX,1);
            
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
            SR{p}.x_all(:,:,j,i) = cohx(:,1:2:end); 
            SR{p}.y_all(:,:,j,i) = cohy(:,2:2:end);
            
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
    SR{p}.x = squeeze(mean(SR{p}.x_all,4));
    SR{p}.y = squeeze(mean(SR{p}.y_all,4));
    SR{p}.xSE = squeeze(std(SR{p}.x_all,[],4)/sqrt(allSubj(p)));
    SR{p}.ySE = squeeze(std(SR{p}.y_all,[],4)/sqrt(allSubj(p)));
    
    RR{p}.x = squeeze(mean(RR{p}.x_all,3));
    RR{p}.y = squeeze(mean(RR{p}.y_all,3));
    RR{p}.xSE = squeeze(std(RR{p}.x_all,[],3)/sqrt(allSubj(p)));
    RR{p}.ySE = squeeze(std(RR{p}.y_all,[],3)/sqrt(allSubj(p)));
    
    x = squeeze(mean(SR{p}.x_all,1));
    y = squeeze(mean(SR{p}.y_all,1));
    
    % total.x_all = 1 - x;
    % total.y_all = 1 - y;
    
    NL{p}.x_all = RR{p}.x_all - x;
    NL{p}.y_all = RR{p}.y_all - y;
    
    NL{p}.x = squeeze(mean(NL{p}.x_all,3));
    NL{p}.y = squeeze(mean(NL{p}.y_all,3));
    NL{p}.xSE = squeeze(std(NL{p}.x_all,[],3)/sqrt(allSubj(p)));
    NL{p}.ySE = squeeze(std(NL{p}.y_all,[],3)/sqrt(allSubj(p)));
    
    % proportion.x_all = NL.x_all./total.x_all;
    % proportion.y_all = NL.y_all./total.y_all;
    %
    % proportion.x = squeeze(mean(proportion.x_all,3));
    % proportion.y = squeeze(mean(proportion.y_all,3));
    % proportion.xSE = squeeze(std(proportion.x_all,[],3)/sqrt(allSubj(p)));
    % proportion.ySE = squeeze(std(proportion.y_all,[],3)/sqrt(allSubj(p)));
end

% generate figures
figure(13); clf
for j = 1:3
    % plot Figure 4B
    subplot(2,3,j); hold on
    for i = 1:Nblock
        for k = 1:Nfreq
            plotIdx = Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials;
            s = shadedErrorBar(plotIdx,SR{j}.x(:,k,i),SR{j}.xSE(:,k,i));
            editErrorBar(s,col(k,:),1.5);
        end
    end
    title(graph_names{j})
    set(gca,'TickDir','out')
    if j == 1
        ylabel('Coherence for x-target freqs')
    end
    xticks(1:5:16)
    xticklabels({'Baseline', 'Early', 'Late', 'Flip'})
    yticks(0:0.25:1)
    axis([1 Ntrials*Nblock 0 1])
    
    % plot Figure 4-supplement 1B
    subplot(2,3,j+3); hold on
    for i = 1:Nblock
        for k = 1:Nfreq
            plotIdx = Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials;
            s = shadedErrorBar(plotIdx,SR{j}.y(:,k,i),SR{j}.ySE(:,k,i));
            editErrorBar(s,col(k,:),1.5);
        end
    end
    set(gca,'TickDir','out')
    if j == 1
        ylabel('Coherence for y-target freqs')
    end
    xticks(1:5:16)
    xticklabels({'Baseline', 'Early', 'Late', 'Flip'})
    yticks(0:0.25:1)
    axis([1 Ntrials*Nblock 0 1])
end

% figure(50); clf
% for j = 1:3
%     subplot(2,3,j); hold on
%     for k = 1:Nfreq
%         s = shadedErrorBar(1:Nblock,NL{j}.x(k,:),NL{j}.xSE(k,:));
%         editErrorBar(s,col(k,:),1.5);
%     end
%     if j == 1
%         title('Rotation')
%     else
%         title('Mirror-Reversal')
%     end
%     set(gca,'TickDir','out')
%     if j == 1
%         ylabel('Coherence for x-target freqs')
%     end
%     xticks([1 2 5 6])
%     xticklabels({'Baseline', 'Early', 'Late', 'Post'})
%     yticks(0:0.25:1)
%     axis([1 6 -0.2 1])
%     
%     subplot(2,3,j+3); hold on
%     for k = 1:Nfreq
%         s = shadedErrorBar(1:Nblock,NL{j}.y(k,:),NL{j}.ySE(k,:));
%         editErrorBar(s,col(k,:),1.5);
%     end
%     set(gca,'TickDir','out')
%     xlabel('Trial Number')
%     if j == 1
%         ylabel('Coherence for y-target freqs')
%     end
%     xticks([1 2 5 6])
%     xticklabels({'Baseline', 'Early', 'Late', 'Post'})
%     yticks(0:0.25:1)
%     axis([1 6 -0.2 1])
% end
end