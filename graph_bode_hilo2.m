function graph_bode_hilo2(data, groups, graph_name, gblocks)

    group_names = {'Rotation','Mirror Reversal'};
    m = figure('Name','High vs. Low Frequency (X Y average)','NumberTitle','off','DefaultTextFontName','Arial','DefaultAxesFontName','Arial'); 
    n = figure('Name','High vs. Low Frequency (X Y separate)','NumberTitle','off','DefaultTextFontName','Arial','DefaultAxesFontName','Arial');
%     o = figure('Name','High vs. Low Frequency (all trajectories)','NumberTitle','off','DefaultTextFontName','Arial','DefaultAxesFontName','Arial');
    p = figure('Name','Aftereffects','NumberTitle','off','DefaultTextFontName','Arial','DefaultAxesFontName','Arial');
    col = [0 0.4470 0.7410];
    col2 = [0.8500 0.3250 0.0980];
    leg = {'Low frequencies','High frequencies'};
    
    for i = 1:length(groups)
        amp_lowx = mean(data.(groups{i}).avg.x_x.d.amplitude(gblocks,1:2),2); %mean amplitude of 2 lowest freqs
        amp_lowy = mean(data.(groups{i}).avg.y_y.d.amplitude(gblocks,1:2),2);
        amp_highx = mean(data.(groups{i}).avg.x_x.d.amplitude(gblocks,6:7),2); %mean amplitude of 2 highest freqs
        amp_highy = mean(data.(groups{i}).avg.y_y.d.amplitude(gblocks,6:7),2);
        
        err_lowx = mean(data.(groups{i}).avg.x_x.d.amp_err(gblocks,1:2,:),2); %mean error of 2 lowest freqs
        err_lowy = mean(data.(groups{i}).avg.y_y.d.amp_err(gblocks,1:2,:),2);
        err_highx = mean(data.(groups{i}).avg.x_x.d.amp_err(gblocks,6:7,:),2); %mean error of 2 highest freqs
        err_highy = mean(data.(groups{i}).avg.y_y.d.amp_err(gblocks,6:7,:),2);
        
        a = [amp_lowx amp_lowy amp_highx amp_highy];
        b = a;
        a = 20*log10(a);
        e2 = [err_lowx err_lowy err_highx err_highy];
        e2 = permute(e2,[3 1 2]);
        e = [e2(2,:,:); e2(1,:,:)];
        x = 1:length(gblocks);

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
        plot(0:1,50:51,'Color',col,'LineWidth',3); hold on;
        plot(0:1,50:51,'Color',col2,'LineWidth',3);
        r = shadedErrorBar(x,mean(a(:,1:2),2),mean(e(:,:,1:2),3),'lineProps','-ob');
        editErrorBar(r,col,3);
        s = shadedErrorBar(x,mean(a(:,3:4),2),mean(e(:,:,3:4),3),'lineProps','-or');
        editErrorBar(s,col2,3);
        patch([1.95 3.05 3.05 1.95],[-80 -80 10 10],'g','FaceAlpha',0.1,'EdgeColor','none');
        patch([4.7 5 5 4.7],[-80 -80 10 10],'c','FaceAlpha',0.1,'EdgeColor','none');
        text(2.5,-2,'ON','HorizontalAlignment','center','FontSize',18);
        text(4.85,-2,'ON','HorizontalAlignment','center','FontSize',18);
        set(gca,'xtick',1:length(x),'xticklabel',graph_name,'ytick',-40:10:0,'FontSize',16,'box','off');
        title(group_names{i},'FontSize',30);
        if i == 1
            ylabel('Magnitude (dB)','FontSize',22)
        end
        ylim([-40 0]);
        legend(leg,'Position',[0.16 0.14 0.1 0.075],'FontSize',16);
        
%% plot x and y separately for high vs low freqs
        figure(n)
        subplot(2,2,i)
        r = shadedErrorBar(x,a(:,1),e(:,:,1),'lineProps','-o'); hold on;
        editErrorBar(r,col,3);
        s = shadedErrorBar(x,a(:,3),e(:,:,3),'lineProps','-o');
        editErrorBar(s,col2,3);
        patch([1.95 3.05 3.05 1.95],[-80 -80 10 10],'g','FaceAlpha',0.1,'EdgeColor','none');
        patch([4.7 5 5 4.7],[-80 -80 10 10],'c','FaceAlpha',0.1,'EdgeColor','none');
        text(2.5,0,'ON','HorizontalAlignment','center','FontSize',18);
        text(4.85,0,'ON','HorizontalAlignment','center','FontSize',18);
        set(gca,'xtick',1:length(gblocks),'xticklabel',graph_name,'ytick',-40:10:0,'FontSize',16,'box','off');
        title(group_names{i},'FontSize',30);
        if i == 1
            ylabel('X Magnitude (dB)','FontSize',22);
        end
        ylim([-50 5]);
        
        subplot(2,2,i+2)
        plot(0:1,50:51,'Color',col,'LineWidth',3); hold on;
        plot(0:1,50:51,'Color',col2,'LineWidth',3);
        r = shadedErrorBar(x,a(:,2),e(:,:,2),'lineProps','-o');
        editErrorBar(r,col,3);
        s = shadedErrorBar(x,a(:,4),e(:,:,4),'lineProps','-o');
        editErrorBar(s,col2,3);
        patch([1.95 3.05 3.05 1.95],[-80 -80 10 10],'g','FaceAlpha',0.1,'EdgeColor','none');
        patch([4.7 5 5 4.7],[-80 -80 10 10],'c','FaceAlpha',0.1,'EdgeColor','none');
        text(2.5,0,'ON','HorizontalAlignment','center','FontSize',18);
        text(4.85,0,'ON','HorizontalAlignment','center','FontSize',18);
        set(gca,'xtick',1:length(gblocks),'xticklabel',graph_name,'ytick',-40:10:0,'FontSize',16,'box','off');
        if i == 1
            ylabel('Y Magnitude (dB)','FontSize',22);
        end
        ylim([-50 5]);
        legend(leg,'FontSize',18,'Position',[0.17 0.12 0.1 0.1]);

%% plot each trial individually
%         a2 = [mean(a2(:,1:2),2) mean(a2(:,3:4),2)];
%         figure(o)
%         subplot(2,1,i)
%         plot(a2);
%         title(group_names{i},'FontSize',15);
%         set(gca,'xtick',1:8:size(a2,1),'xticklabel',graph_name,'FontSize',12);
%         ylabel('Magnitude (dB)','FontSize',12);
%         ylim([-40 5]);
%         legend(leg,'Position',[0 0.4 0.1 0.075],'FontSize',12);
%         grid on;

%% aftereffects
        aftereffect = b(1,:) - b(4,:);
        total =  b(1,:) - b(2,:);
        leg2 = {'Low X','Low Y','High X','High Y'};
        
        figure(p)
        subplot(1,2,i)
        bar(100*(aftereffect./total));
        set(gca,'FontSize',22);
        title(group_names{i},'FontSize',30)
        ylabel('Aftereffect (%)');
        xticklabels(leg2);
        ylim([0 100]);
    end
end

function output = editErrorBar(x, col, width)
    x.mainLine.Color = col;
    x.mainLine.LineWidth = width;
    x.patch.FaceColor = col;
    x.patch.EdgeColor = col+(1-col)*0.55;
    set(x.edge,'Color',col)
    
    output = x;
end