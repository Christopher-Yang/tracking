function graph_lag(data,gblocks,graph_name,output)

    group_names = {'Rotation','Mirror Reversal'};
    dat = data{end}.(output);
    f_x = dat.x_x.freqs;
    f_y = dat.y_y.freqs;
    col = lines;
    col = col(1:7,:);
    Nblock = size(dat.x_x.fft,1);
    
    delay.x_x = -1000*dat.x_x.phase./(repmat(f_x,[Nblock,1])*2*pi);
    delay.y_y = -1000*dat.y_y.phase./(repmat(f_y,[Nblock,1])*2*pi);
    delay.x_y = -1000*dat.x_y.phase./(repmat(f_x,[Nblock,1])*2*pi);
    delay.y_x = -1000*dat.y_x.phase./(repmat(f_y,[Nblock,1])*2*pi);
    delay.tc = -1000*mean(cat(3,delay.x_x(gblocks,:),delay.y_y(gblocks,:)),3);
    delay.hc = -1000*mean(cat(3,delay.x_y(gblocks,:),delay.y_x(gblocks,:)),3);
    err.x_x = permute(1000*dat.x_x.phase_err./(repmat(f_x,[Nblock,1,2])*360),[3 2 1]);
    err.y_y = permute(1000*dat.y_y.phase_err./(repmat(f_y,[Nblock,1,2])*360),[3 2 1]);
    err.x_y = permute(1000*dat.x_y.phase_err./(repmat(f_x,[Nblock,1,2])*360),[3 2 1]);
    err.y_x = permute(1000*dat.y_x.phase_err./(repmat(f_y,[Nblock,1,2])*360),[3 2 1]);
            
    figure
    subplot(2,1,1); hold on
    plot(f_x,delay.x_x(gblocks(1),:)','--o','Color',col(1,:),'LineWidth',1,'MarkerFaceColor',col(1,:))
    plot(f_x,delay.x_x(gblocks(2),:)','-o','Color',col(1,:),'LineWidth',1,'MarkerFaceColor',col(1,:))
    plot(f_x,delay.x_x(gblocks(3),:)','--o','Color',col(2,:),'LineWidth',1,'MarkerFaceColor',col(2,:))
    plot(f_x,delay.x_x(gblocks(4),:)','-o','Color',col(2,:),'LineWidth',1,'MarkerFaceColor',col(2,:))
%     s = shadedErrorBar(f_x,delay.x_x(gblocks(1),:)',err.x_x(:,:,gblocks(1)),'lineProps','--o'); hold on
%     editErrorBar(s,col(1,:),2);
%     s = shadedErrorBar(f_x,delay.x_x(gblocks(2),:)',err.x_x(:,:,gblocks(2)),'lineProps','-o'); hold on
%     editErrorBar(s,col(1,:),2);
%     s = shadedErrorBar(f_x,delay.x_x(gblocks(3),:)',err.x_x(:,:,gblocks(3)),'lineProps','--o'); hold on
%     editErrorBar(s,col(2,:),2);
%     s = shadedErrorBar(f_x,delay.x_x(gblocks(4),:)',err.x_x(:,:,gblocks(4)),'lineProps','-o'); hold on
%     editErrorBar(s,col(2,:),2);
    set(gca,'Xscale','log','Xgrid','on','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','TickDir','out','LineWidth',1)
    title('X frequencies')
    ylabel('Response Lag (ms)')
    axis([0.09 2.3 0 600])

    subplot(2,1,2); hold on
    plot(f_y,delay.y_y(gblocks(1),:)','--o','Color',col(1,:),'LineWidth',1,'MarkerFaceColor',col(1,:))
    plot(f_y,delay.y_y(gblocks(2),:)','-o','Color',col(1,:),'LineWidth',1,'MarkerFaceColor',col(1,:))
    plot(f_y,delay.y_y(gblocks(3),:)','--o','Color',col(2,:),'LineWidth',1,'MarkerFaceColor',col(2,:))
    plot(f_y,delay.y_y(gblocks(4),:)','-o','Color',col(2,:),'LineWidth',1,'MarkerFaceColor',col(2,:))
%     s = shadedErrorBar(f_y,delay.y_y(gblocks(1),:)',err.y_y(:,:,gblocks(1)),'lineProps','--o'); hold on
%     editErrorBar(s,col(1,:),2);
%     s = shadedErrorBar(f_y,delay.y_y(gblocks(2),:)',err.y_y(:,:,gblocks(2)),'lineProps','-o'); hold on
%     editErrorBar(s,col(1,:),2);
%     s = shadedErrorBar(f_y,delay.y_y(gblocks(3),:)',err.y_y(:,:,gblocks(3)),'lineProps','--o'); hold on
%     editErrorBar(s,col(2,:),2);
%     s = shadedErrorBar(f_y,delay.y_y(gblocks(4),:)',err.y_y(:,:,gblocks(4)),'lineProps','-o'); hold on
%     editErrorBar(s,col(2,:),2);
    set(gca,'Xscale','log','Xgrid','on','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','TickDir','out','LineWidth',1)
    ylabel('Response Lag (ms)')
    xlabel('Frequency (Hz)')
    axis([0.09 2.3 0 600])
    legend(graph_name{gblocks},'Position',[0.78 0.77 0.1 0.1],'FontSize',12)
end