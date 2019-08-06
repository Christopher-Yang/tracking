function graph_MSE(data,block_name,graph_name)

    col = lines;
    col = [col(1:7,:) repelem(0.1,7)'];
    col2 = [160 82 45
            46 139 87
            65 105 225]./255;
    
    Nsubj = length(data)-1;
    Nblocks = length(block_name);
    Ntrials = length(data{1}.(block_name{1}).MSE);
    MSE_subj = NaN(Nblocks, Nsubj);
    MSE_se = NaN(Nblocks, Nsubj);
    MSE_full = NaN(Nblocks*Ntrials, Nsubj);
    for j = 1:Nsubj
        x = NaN(Nblocks*Ntrials, 1);
        for k = 1:Nblocks
            MSE_subj(k,j) = mean(data{j}.(block_name{k}).MSE);
            MSE_se(k,j) = std(data{j}.(block_name{k}).MSE)/sqrt(length(data{j}.(block_name{k}).MSE));
            x((k-1)*Ntrials+1:(k-1)*Ntrials+Ntrials) = data{j}.(block_name{k}).MSE;
        end
        MSE_full(:,j) = x;
    end
    MSE.se2 = std(MSE_subj,0,2)/sqrt(size(MSE_subj,2)); % se2 is the standard error of average across subjects
    MSE.subjects = MSE_subj;
    MSE.se = MSE_se; % MSE_se is the standard error across trials within individuals
    MSE.full = MSE_full;
    MSE.full_se = std(MSE_full,0,2)/sqrt(size(MSE_full,2));
%     figure('Name', 'Mean Squared Error (by condition)','NumberTitle','off')
%     errorbar(mean(MSE.rot.subjects,2), MSE.rot.se2,'r','LineWidth',1.5); hold on; 
%     errorbar(mean(MSE.rot_i.subjects,2), MSE.rot_i.se2, 'b','LineWidth',1.5);
%     set(gca,'xtick',1:6,'xticklabel',graph_name);
%     ylabel('Mean Squared Error (m^2)'); 
%     legend('90 deg','mirror reversal');
%     grid on;
    
    figure(1); clf; hold on
    errorbar(mean(MSE.subjects,2), 2*MSE.se2,'r','LineWidth',1);
    plot(MSE.subjects(:,1),'Color',[1 0 0 0.25],'LineWidth',0.5) 
    plot(MSE.subjects(:,2:end),'Color',[1 0 0 0.25],'LineWidth',0.5)
    set(gca,'xtick',1:Nblocks,'ytick',0:5:10,'TickDir','out')
    ylabel('Mean Squared Error (m^2)') 
    box off

    figure(2); clf; hold on
    for i = 1:length(block_name)
        s = shadedErrorBar(Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials,mean(MSE.full(Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials,:),2),2*MSE.full_se(Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials));
        editErrorBar(s,col2(1,:),1);
        plot(Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials,MSE.full(Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials,:),'Color',[col2(3,:) 0.25],'LineWidth',0.5)
    end
    set(gca,'xtick',1:Ntrials:length(MSE.full),'ytick',0:0.003:0.012,'TickDir','out')
%     axis([1 Ntrials*Nblocks 0 0.012])
    xlabel('Tracking Cycle Number')
    ylabel('Mean Squared Error (m^2)'); 
    box off
end