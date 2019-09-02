function graph_lag(data,gblocks,graph_name,output)

    dat = data{end}.(output);
    f_x = dat.x_x.freqs;
    f_y = dat.y_y.freqs;
    col = lines;
    col = col(1:7,:);
    Nblock = size(dat.x_x.fft,1);
    Nsubj = length(data)-1;
    
    names = {'x_x','y_y','x_y','y_x'};
    
    delay.x_x = -1000*dat.x_x.phase./(repmat(f_x,[Nblock,1])*2*pi);
    delay.y_y = -1000*dat.y_y.phase./(repmat(f_y,[Nblock,1])*2*pi);
    delay.x_y = -1000*dat.x_y.phase./(repmat(f_x,[Nblock,1])*2*pi);
    delay.y_x = -1000*dat.y_x.phase./(repmat(f_y,[Nblock,1])*2*pi);
    err.x_x = permute(1000*dat.x_x.phase_err./(repmat(f_x,[Nblock,1,2])*360),[3 2 1]);
    err.y_y = permute(1000*dat.y_y.phase_err./(repmat(f_y,[Nblock,1,2])*360),[3 2 1]);
    err.x_y = permute(1000*dat.x_y.phase_err./(repmat(f_x,[Nblock,1,2])*360),[3 2 1]);
    err.y_x = permute(1000*dat.y_x.phase_err./(repmat(f_y,[Nblock,1,2])*360),[3 2 1]);
    
    delay_all.x_x = permute(-1000*dat.x_x.all_phase./(repmat(f_x,[Nblock,1,Nsubj])*360),[2 3 1]);
    delay_all.y_y = permute(-1000*dat.y_y.all_phase./(repmat(f_y,[Nblock,1,Nsubj])*360),[2 3 1]);
    delay_all.x_y = permute(-1000*dat.x_y.all_phase./(repmat(f_x,[Nblock,1,Nsubj])*360),[2 3 1]);
    delay_all.y_x = permute(-1000*dat.y_x.all_phase./(repmat(f_y,[Nblock,1,Nsubj])*360),[2 3 1]);
    
    figure
    subplot(2,1,1); hold on
    for i = 1:length(gblocks)
        if i == 1 || i == 2
%             s = shadedErrorBar(f_x,delay.x_x(gblocks(i),:)',err.x_x(:,:,gblocks(i)),'lineProps','-o');
%             editErrorBar(s,col(i,:),2);
            plot(f_x,delay.x_x(gblocks(i),:)','-o','Color',col(i,:),'linewidth',1,'MarkerFaceColor',col(i,:))
            plot(f_x,delay_all.x_x(:,:,i),'Color',[col(i,:) 0.4])
        else
%             s = shadedErrorBar(f_x,delay.x_y(gblocks(i),:)',err.x_y(:,:,gblocks(i)),'lineProps','-o');
%             editErrorBar(s,col(i,:),2);
            plot(f_y,delay.y_x(gblocks(i),:)','-o','Color',col(i,:),'linewidth',1,'MarkerFaceColor',col(i,:))
            plot(f_y,delay_all.y_x(:,:,i),'Color',[col(i,:) 0.4])
        end
    end
    set(gca,'box','off','TickDir','out')
    ylabel(['X ',output,' response lag (ms)'])
    xlim([0.09 2.3])

    subplot(2,1,2); hold on
    for i = 1:length(gblocks)
        plot([-5 -6],[0 0],'-o','linewidth',1,'MarkerFaceColor',col(i,:))
    end
    for i = 1:length(gblocks)
        if i == 1 || i == 2
%             s = shadedErrorBar(f_y,delay.y_y(gblocks(i),:)',err.y_y(:,:,gblocks(i)),'lineProps','-o');
%             editErrorBar(s,col(i,:),2);
            plot(f_y,delay.y_y(gblocks(i),:)','-o','Color',col(i,:),'linewidth',1,'MarkerFaceColor',col(i,:))
            plot(f_y,delay_all.y_y(:,:,i),'Color',[col(i,:) 0.4])
        else
%             s = shadedErrorBar(f_y,delay.y_x(gblocks(i),:)',err.y_x(:,:,gblocks(i)),'lineProps','-o');
%             editErrorBar(s,col(i,:),2);
            plot(f_x,delay.x_y(gblocks(i),:)','-o','Color',col(i,:),'linewidth',1,'MarkerFaceColor',col(i,:))
            plot(f_x,delay_all.x_y(:,:,i),'Color',[col(i,:) 0.4])
        end
    end
    set(gca,'box','off','TickDir','out')
    ylabel(['Y ',output,' response lag (ms)'])
    xlabel('Frequency (Hz)')
    xlim([0.09 2.3])
    legend(graph_name{gblocks})
end