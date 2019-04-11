function graph_bode_hilo(data, groups, graph_name, gblocks)

    rng(1);
    col = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330; 0.6350 0.0780 0.1840];
    group_names = {'Rotation','Mirror Reversal'};
    m = figure('Name','High vs. Low Frequency (X Y average)','NumberTitle','off','DefaultTextFontName','Arial','DefaultAxesFontName','Arial'); 
%     n = figure('Name','High vs. Low Frequency (X Y separate)','NumberTitle','off','DefaultTextFontName','Arial','DefaultAxesFontName','Arial'); 
%     o = figure('Name','High vs. Low Frequency (all trajectories)','NumberTitle','off','DefaultTextFontName','Arial','DefaultAxesFontName','Arial');
%     p = figure('Name','Aftereffects','NumberTitle','off','DefaultTextFontName','Arial','DefaultAxesFontName','Arial');
    leg = {'Low frequencies','High frequencies'};
%     names = {'x_x','y_y','x_y','y_x'};
    gains = [0.05 0.1 0.25 0.5 1];
    Yt = 20*log10(gains);
    Ytlab = num2cell(gains);

    for i = 1:length(groups)
%         all = struct('x_x',NaN(6,2,10),'y_y',NaN(6,2,10),'x_y',NaN(6,2,10),'y_x',NaN(6,2,10));
        all = NaN(6,2,10,2);
%         avg = struct('x_x',NaN(6,2),'y_y',NaN(6,2),'x_y',NaN(6,2),'y_x',NaN(6,2));
%         boot = struct('x_x',NaN(6,2,1000),'y_y',NaN(6,2,1000),'x_y',NaN(6,2,1000),'y_x',NaN(6,2,1000));
        boot = NaN(6,2,1000);
%         for j = 1:length(names)
            all(:,1,:,1) = mean(data.(groups{i}).avg.x_x.d.all_amp(:,1:2,:),2);
            all(:,2,:,1) = mean(data.(groups{i}).avg.x_x.d.all_amp(:,6:7,:),2);
            all(:,1,:,2) = mean(data.(groups{i}).avg.y_y.d.all_amp(:,1:2,:),2);
            all(:,2,:,2) = mean(data.(groups{i}).avg.y_y.d.all_amp(:,6:7,:),2);
%             all = 20*log10(mean(all,4));
            all = mean(all,4);
            avg = 20*log10(mean(all,3));
            for k = 1:1000
                z = datasample(all,10,3);
                boot(:,:,k) = mean(z,3);
            end
            boot = 20*log10(sort(boot,3));
%         end

%         amp_lowx = mean(data.(groups{i}).avg.x_x.d.amplitude(gblocks,1:2),2); %mean amplitude of 2 lowest freqs
%         amp_lowy = mean(data.(groups{i}).avg.y_y.d.amplitude(gblocks,1:2),2);
%         amp_highx = mean(data.(groups{i}).avg.x_x.d.amplitude(gblocks,6:7),2); %mean amplitude of 2 highest freqs
%         amp_highy = mean(data.(groups{i}).avg.y_y.d.amplitude(gblocks,6:7),2);
        
        err_lowfreq1 = boot(:,1,end-25) - avg(:,1);
        err_lowfreq2 = avg(:,1) - boot(:,1,26);
        err_highfreq1 = boot(:,2,end-25) - avg(:,2);
        err_highfreq2 = avg(:,2) - boot(:,2,26);

%         err_lowx = mean(data.(groups{i}).avg.x_x.d.amp_err(gblocks,1:2,:),2); %mean error of 2 lowest freqs
%         err_lowy = mean(data.(groups{i}).avg.y_y.d.amp_err(gblocks,1:2,:),2);
%         err_highx = mean(data.(groups{i}).avg.x_x.d.amp_err(gblocks,6:7,:),2); %mean error of 2 highest freqs
%         err_highy = mean(data.(groups{i}).avg.y_y.d.amp_err(gblocks,6:7,:),2);
        
        a = avg';
        b = a;
        e2 = [err_lowfreq1 err_lowfreq2 err_highfreq1 err_highfreq2];
%         e2 = permute(e2,[3 1 2]);
        e = e2';
%         e = [e2(2,:,:); e2(1,:,:)];
        x = 1:length(gblocks);
%          x = 1:4;
        
        amp_lowx2 = squeeze(mean(data.(groups{i}).avg.x_x.d2.amplitude(:,1:2,:),2))';
        amp_lowy2 = squeeze(mean(data.(groups{i}).avg.y_y.d2.amplitude(:,1:2,:),2))';
        amp_highx2 = squeeze(mean(data.(groups{i}).avg.x_x.d2.amplitude(:,6:7,:),2))';
        amp_highy2 = squeeze(mean(data.(groups{i}).avg.y_y.d2.amplitude(:,6:7,:),2))';
        
        amp_lowx2 = amp_lowx2(:);
        amp_lowy2 = amp_lowy2(:);
        amp_highx2 = amp_highx2(:);
        amp_highy2 = amp_highy2(:);
        
        a2 = [amp_lowx2 amp_lowy2 amp_highx2 amp_highy2];
        a2 = 20*log10(a2);
        
%% plot x and y together for high vs low freqs
        figure(m)
        subplot(1,2,i)
%         plot(0:1,50:51,'Color',col,'LineWidth',3); hold on;
%         plot(0:1,50:51,'Color',col2,'LineWidth',3);
        
        r = shadedErrorBar(x,a(1,:),e(1:2,:),'lineProps','-ob');
        editErrorBar(r,[158 0 93]/255,3);
        s = shadedErrorBar(x,a(2,:),e(3:4,:),'lineProps','-or');
        editErrorBar(s,[96 56 19]/255,3);
        
        patch([1.95 5.05 5.05 1.95],[-45 -45 10 10],'g','FaceAlpha',0.1,'EdgeColor','none');
        text(3.5,-2,'PERTURBATION ON','HorizontalAlignment','center','FontSize',18);
        set(gca,'xtick',[1:2 5:6],'xticklabel',graph_name,'ytick',Yt,'yticklabel',Ytlab,'FontSize',12,'TickDir','out','LineWidth',1,'box','off');
        title(group_names{i},'FontSize',30);
        if i == 1
            ylabel('Target to cursor gain (cm/cm)','FontSize',22)
        end
        ylim([-30 0]);
        legend(leg,'Position',[0.17 0.22 0.1 0.075],'FontSize',18);
        pbaspect([1 1 1]);
        
%% plot x and y separately for high vs low freqs
%         figure(n)
%         subplot(2,2,i)
%         r = shadedErrorBar(x,a(:,1),e(:,:,1),'lineProps','-o'); hold on;
%         editErrorBar(r,col(6,:),3);
%         s = shadedErrorBar(x,a(:,3),e(:,:,3),'lineProps','-o');
%         editErrorBar(s,col(7,:),3);
%         patch([1.95 5.05 5.05 1.95],[-45 -45 10 10],'g','FaceAlpha',0.1,'EdgeColor','none');
%         text(3.5,0,'ON','HorizontalAlignment','center','FontSize',18);
%         set(gca,'xtick',[1:2 5:6],'xticklabel',graph_name,'ytick',-40:10:0,'FontSize',12,'TickDir','out','LineWidth',1,'box','off');
%         title(group_names{i},'FontSize',30);
%         if i == 1
%             ylabel('X Magnitude (dB)','FontSize',22);
%         end
%         ylim([-40 3]);
%         
%         subplot(2,2,i+2)
%         plot(0:1,50:51,'Color',col(6,:),'LineWidth',3); hold on;
%         plot(0:1,50:51,'Color',col(7,:),'LineWidth',3);
%         r = shadedErrorBar(x,a(:,2),e(:,:,2),'lineProps','-o');
%         editErrorBar(r,col(6,:),3);
%         s = shadedErrorBar(x,a(:,4),e(:,:,4),'lineProps','-o');
%         editErrorBar(s,col(7,:),3);
%         patch([1.95 5.05 5.05 1.95],[-45 -45 10 10],'g','FaceAlpha',0.1,'EdgeColor','none');
%         text(3.5,0,'ON','HorizontalAlignment','center','FontSize',18);
%         set(gca,'xtick',[1:2 5:6],'xticklabel',graph_name,'ytick',-40:10:0,'FontSize',16,'TickDir','out','LineWidth',1,'box','off');
%         if i == 1
%             ylabel('Y Magnitude (dB)','FontSize',22);
%         end
%         ylim([-40 3]);
%         legend(leg,'FontSize',18,'Position',[0.17 0.12 0.1 0.1]);

%% individual trials
%          a2 = [mean(a2(:,1:2),2) mean(a2(:,3:4),2)];
%          figure(n)
%          subplot(2,1,i)
%          plot(a2);
%          title(group_names{i},'FontSize',15);
%          set(gca,'xtick',1:8:size(a2,1),'xticklabel',graph_name,'FontSize',12);
%          ylabel('Magnitude (dB)','FontSize',12);
%          ylim([-40 5]);
%          legend(leg,'Position',[0 0.4 0.1 0.075],'FontSize',12);
%          grid on;

% %% aftereffects
% %         aftereffect = b(1,:) - b(6,:);
% %         total =  b(1,:) - b(2,:);
% 
%         c = [mean(b(:,1:2),2) mean(b(:,3:4),2)];
%         aftereffect = c(1,:) - c(6,:);
%         total = c(1,:) - c(2,:);
%         
% %         leg2 = {'Low X','Low Y','High X','High Y'};
%         leg2 = {'Low Frequencies','High Frequencies'};
% 
%         figure(p)
%         subplot(1,2,i)
%         bar(100*(aftereffect./total));
%         set(gca,'FontSize',22)
%         title(group_names{i},'FontSize',30)
%         ylabel('Aftereffect (%)');
%         xticklabels(leg2);
%         ylim([0 100]);
    end
end