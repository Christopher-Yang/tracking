function graph_lag(data,gblocks,graph_name)

    group_names = {'Rotation','Mirror Reversal'};
    f_x = data.(groups{1}).avg.x_x.freqs;
    f_y = data.(groups{1}).avg.y_y.freqs;
    col = lines;
    col = col(1:7,:);
    Nblock = size(data.(groups{1}).avg.x_x.fft,1);
    
    for i = 1:length(groups)
        delay{i}.x_x = -1000*data.(groups{i}).avg.x_x.phase./(repmat(f_x,[Nblock,1])*2*pi);
        delay{i}.y_y = -1000*data.(groups{i}).avg.y_y.phase./(repmat(f_y,[Nblock,1])*2*pi);
        delay{i}.x_y = -1000*data.(groups{i}).avg.x_y.phase./(repmat(f_x,[Nblock,1])*2*pi);
        delay{i}.y_x = -1000*data.(groups{i}).avg.y_x.phase./(repmat(f_y,[Nblock,1])*2*pi);
        delay{i}.tc = -1000*mean(cat(3,delay{i}.x_x(gblocks,:),delay{i}.y_y(gblocks,:)),3);
        delay{i}.hc = -1000*mean(cat(3,delay{i}.x_y(gblocks,:),delay{i}.y_x(gblocks,:)),3);
        err{i}.x_x = permute(1000*data.(groups{i}).avg.x_x.phase_err./(repmat(f_x,[Nblock,1,2])*360),[3 2 1]);
        err{i}.y_y = permute(1000*data.(groups{i}).avg.y_y.phase_err./(repmat(f_y,[Nblock,1,2])*360),[3 2 1]);
        err{i}.x_y = permute(1000*data.(groups{i}).avg.x_y.phase_err./(repmat(f_x,[Nblock,1,2])*360),[3 2 1]);
        err{i}.y_x = permute(1000*data.(groups{i}).avg.y_x.phase_err./(repmat(f_y,[Nblock,1,2])*360),[3 2 1]);
    end 
    
%     o = figure('Name','Response Lag','NumberTitle','off','DefaultTextFontName','Arial','DefaultAxesFontName','Arial');
%     p = figure('Name','Response Lag (simple)','NumberTitle','off','DefaultTextFontName','Arial','DefaultAxesFontName','Arial');
    q = figure('Name','Response Lag (congruent)','NumberTitle','off','DefaultTextFontName','Arial','DefaultAxesFontName','Arial');
    
    for i = 1:length(groups)
%         figure(o)
%         subplot(2,2,(i-1)*2 + 1)
%         semilogx(f_x,delay{i}.x_x','-o','LineWidth',1.5)
%         set(gca,'FontSize',12,'TickDir','out','LineWidth',1)
%         title([group_names{i},' (X)'],'FontSize',15)
%         ylabel('Response Lag (ms)','FontSize',15)
%         xlabel('Frequency (Hz)','FontSize',15)
%         axis([0.09 2.3 300 1500])
%         grid on
%         
%         subplot(2,2,(i-1)*2 + 2)
%         semilogx(f_y,delay{i}.y_y','-o','LineWidth',1.5)
%         set(gca,'FontSize',12,'TickDir','out','LineWidth',1)
%         title([group_names{i},' (Y)'],'FontSize',15)
%         ylabel('Response Lag (ms)','FontSize',15)
%         xlabel('Frequency (Hz)','FontSize',15)
%         axis([0.09 2.3 300 1500]);
%         grid on
%         legend(graph_name{gblocks},'Position',[0.9 0.45 0.1 0.1],'FontSize',12)
        
%         figure(p)        
%         subplot(2,2,i)
%         for j = 1:length(gblocks)
%             s = shadedErrorBar(f_x,delay{i}.x_x(gblocks(j),:)',err{i}.x_x(:,:,gblocks(j)),'lineProps','-o'); hold on;
%             editErrorBar(s,col(j,:),3);
%         end
%         set(gca,'Xscale','log','FontSize',18,'Xgrid','on','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','TickDir','out','LineWidth',1,'Ytick',300:300:1500)
%         title(group_names{i},'FontSize',30)
%         if i == 1
%             ylabel('X Response Lag (ms)','FontSize',22)
%         end
%         axis([0.09 2.3 300 1500])
%         
%         subplot(2,2,i+2)
%         for j = 1:length(gblocks)
%             s = shadedErrorBar(f_y,delay{i}.y_y(gblocks(j),:)',err{i}.y_y(:,:,gblocks(j)),'lineProps','-o');
%             editErrorBar(s,col(j,:),3);
%         end
%         set(gca,'Xscale','log','FontSize',18,'Xgrid','on','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','TickDir','out','LineWidth',1,'Ytick',300:300:1500)
%         if i == 1
%             ylabel('Y Response Lag (ms)','FontSize',22)
%         end
%         xlabel('Frequency (Hz)','FontSize',22)
%         axis([0.09 2.3 300 1500])
%         legend(graph_name{gblocks},'Position',[0.34 0.8 0.1 0.1],'FontSize',18)
        
        figure(q)
        subplot(2,2,i)
        s = shadedErrorBar(f_x,delay{i}.x_x(gblocks(1),:)',err{i}.x_x(:,:,gblocks(1)),'lineProps','-o'); hold on
        editErrorBar(s,col(1,:),2);
        s = shadedErrorBar(f_y,delay{i}.y_y(gblocks(2),:)',err{i}.y_y(:,:,gblocks(2)),'lineProps','-o');
        editErrorBar(s,col(2,:),2);
        s = shadedErrorBar(f_y,delay{i}.y_y(gblocks(3),:)',err{i}.y_y(:,:,gblocks(3)),'lineProps','-o');
        editErrorBar(s,col(3,:),2);
        s = shadedErrorBar(f_x,delay{i}.x_x(gblocks(4),:)',err{i}.x_x(:,:,gblocks(4)),'lineProps','-o');
        editErrorBar(s,col(4,:),2);
        set(gca,'Xscale','log','FontSize',12,'Xgrid','on','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','TickDir','out','LineWidth',1,'Ytick',300:300:1500)
        title(group_names{i},'FontSize',15)
        ylabel('Response Lag (ms)','FontSize',15)
        axis([0.09 2.3 300 1500])
        
        subplot(2,2,i+2)
        s = shadedErrorBar(f_y,delay{i}.y_y(gblocks(1),:)',err{i}.y_y(:,:,gblocks(1)),'lineProps','-o'); hold on
        editErrorBar(s,col(1,:),2);
        s = shadedErrorBar(f_x,delay{i}.x_x(gblocks(2),:)',err{i}.x_x(:,:,gblocks(2)),'lineProps','-o');
        editErrorBar(s,col(2,:),2);
        s = shadedErrorBar(f_x,delay{i}.x_x(gblocks(3),:)',err{i}.x_x(:,:,gblocks(3)),'lineProps','-o');
        editErrorBar(s,col(3,:),2);
        s = shadedErrorBar(f_y,delay{i}.y_y(gblocks(4),:)',err{i}.y_y(:,:,gblocks(4)),'lineProps','-o');
        editErrorBar(s,col(4,:),2);
        set(gca,'Xscale','log','FontSize',12,'Xgrid','on','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','TickDir','out','LineWidth',1,'Ytick',300:300:1500)
        ylabel('Response Lag (ms)','FontSize',15)
        xlabel('Frequency (Hz)','FontSize',15)
        axis([0.09 2.3 300 1500])
        legend(graph_name{gblocks},'Position',[0.78 0.77 0.1 0.1],'FontSize',12)
    end
end