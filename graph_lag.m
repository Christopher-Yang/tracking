function graph_lag(data,gblocks,graph_name,output)

    dat = data{end}.(output);
    f_x = dat.x_x.freqs;
    f_y = dat.y_y.freqs;
    Nblock = size(dat.x_x.fft,1);
    Nfreq = length(f_x);
    col = copper;
    n = round(64.*((1:Nfreq-1)/(Nfreq-1)));
    col = col([1 n],:);
    
    delay.x_x = -1000*dat.x_x.phase./(repmat(f_x,[Nblock,1])*2*pi);
    delay.y_y = -1000*dat.y_y.phase./(repmat(f_y,[Nblock,1])*2*pi);
    delay.x_y = -1000*dat.x_y.phase./(repmat(f_x,[Nblock,1])*2*pi);
    delay.y_x = -1000*dat.y_x.phase./(repmat(f_y,[Nblock,1])*2*pi);
    err.x_x = permute(1000*dat.x_x.phase_err./(repmat(f_x,[Nblock,1,2])*360),[3 2 1]);
    err.y_y = permute(1000*dat.y_y.phase_err./(repmat(f_y,[Nblock,1,2])*360),[3 2 1]);
    err.x_y = permute(1000*dat.x_y.phase_err./(repmat(f_x,[Nblock,1,2])*360),[3 2 1]);
    err.y_x = permute(1000*dat.y_x.phase_err./(repmat(f_y,[Nblock,1,2])*360),[3 2 1]);
        
%     figure
%     subplot(2,1,1); hold on
%     for i = 1:length(gblocks)
%         if i == 1
% %             s = shadedErrorBar(f_x,delay.x_x(gblocks(i),:)',err.x_x(:,:,gblocks(i)),'lineProps','-o');
% %             editErrorBar(s,col(i,:),2);
%             plot(f_x,delay.x_x(gblocks(i),:)','-o','Color',col(i,:),'linewidth',1);
%         else
% %             s = shadedErrorBar(f_x,delay.x_y(gblocks(i),:)',err.x_y(:,:,gblocks(i)),'lineProps','-o');
% %             editErrorBar(s,col(i,:),2);
%             plot(f_y,delay.y_x(gblocks(i),:)','-o','Color',col(i,:),'linewidth',1);
%         end
%     end
%     set(gca,'Xscale','log','FontSize',12,'Xgrid','on','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','TickDir','out','LineWidth',1)
%     title(group_names,'FontSize',15)
%     ylabel('Response Lag (ms)','FontSize',15)
%     xlim([0.09 2.3])
% 
%     subplot(2,1,2); hold on
%     for i = 1:length(gblocks)
%         if i == 1
% %             s = shadedErrorBar(f_y,delay.y_y(gblocks(i),:)',err.y_y(:,:,gblocks(i)),'lineProps','-o');
% %             editErrorBar(s,col(i,:),2);
%             plot(f_y,delay.y_y(gblocks(i),:)','-o','Color',col(i,:),'linewidth',1);
%         else
% %             s = shadedErrorBar(f_y,delay.y_x(gblocks(i),:)',err.y_x(:,:,gblocks(i)),'lineProps','-o');
% %             editErrorBar(s,col(i,:),2);
%             plot(f_x,delay.x_y(gblocks(i),:)','-o','Color',col(i,:),'linewidth',1);
%         end
%     end
%     set(gca,'Xscale','log','FontSize',12,'Xgrid','on','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','TickDir','out','LineWidth',1)
%     ylabel('Response Lag (ms)','FontSize',15)
%     xlabel('Frequency (Hz)','FontSize',15)
%     xlim([0.09 2.3])
%     legend(graph_name{gblocks},'Position',[0.78 0.77 0.1 0.1],'FontSize',12)
%     
    for i = 1:Nfreq
        leg{i} = ['Freq ',num2str(i)];
    end
    tick = find(contains(graph_name,'('));
    
    figure
    subplot(2,1,1); hold on
    for i = 1:Nfreq
        plot(1:Nblock,delay.x_x(:,i)','-o','Color',col(i,:),'LineWidth',1,'MarkerFaceColor',col(i,:),'MarkerEdgeColor','none')
    end
    set(gca,'TickDir','out')
    xlabel('Block')
    ylabel([output,' X Lag (ms)'])
    xticks(tick)
    xticklabels(graph_name(tick))
    legend(leg)
    
    subplot(2,1,2); hold on
    for i = 1:Nfreq
        plot(1:Nblock,delay.y_y(:,i)','-o','Color',col(i,:),'LineWidth',1,'MarkerFaceColor',col(i,:),'MarkerEdgeColor','none')
    end
    set(gca,'TickDir','out')
    xlabel('Block')
    ylabel([output,' Y Lag (ms)'])
    xticks(tick)
    xticklabels(graph_name(tick))
end