function graph_lag(data,gblocks,graph_name,output)

    group_names = {'Rotation','Mirror Reversal'};
    dat = data{end}.(output);
    f_x = data{end}.freqX;
    f_y = data{end}.freqY;
    Nblock = length(gblocks);
    Nsubj = length(data)-1;
    c = copper;
    col1 = c(1,:);
    col2 = c(end,:);
    col = [linspace(col1(1),col2(1),Nblock)', linspace(col1(2),col2(2),Nblock)', linspace(col1(3),col2(3),Nblock)'];
    
    delay.x_x = -1000*dat.x_x.phase./(repmat(f_x,[Nblock 1])*2*pi);
    delay.y_y = -1000*dat.y_y.phase./(repmat(f_y,[Nblock 1])*2*pi);
    delay.x_y = -1000*dat.x_y.phase./(repmat(f_x,[Nblock 1])*2*pi);
    delay.y_x = -1000*dat.y_x.phase./(repmat(f_y,[Nblock 1])*2*pi);
    delay.x_xAll = permute(-1000*dat.x_x.all_phase./(repmat(f_x,[Nblock 1 Nsubj])*2*pi),[2 1 3]);
    delay.y_yAll = permute(-1000*dat.y_y.all_phase./(repmat(f_y,[Nblock 1 Nsubj])*2*pi),[2 1 3]);
    delay.x_yAll = permute(-1000*dat.x_y.all_phase./(repmat(f_x,[Nblock 1 Nsubj])*2*pi),[2 1 3]);
    delay.y_xAll = permute(-1000*dat.y_x.all_phase./(repmat(f_y,[Nblock 1 Nsubj])*2*pi),[2 1 3]);

%     err.x_x = permute(1000*dat.x_x.phase_err./(repmat(f_x,[Nblock,1,2])*360),[3 2 1]);
%     err.y_y = permute(1000*dat.y_y.phase_err./(repmat(f_y,[Nblock,1,2])*360),[3 2 1]);
%     err.x_y = permute(1000*dat.x_y.phase_err./(repmat(f_x,[Nblock,1,2])*360),[3 2 1]);
%     err.y_x = permute(1000*dat.y_x.phase_err./(repmat(f_y,[Nblock,1,2])*360),[3 2 1]);
    figure
    for j = 1:Nsubj
        subplot(2,Nsubj,j); hold on
%         plot(f_x,delay.x_x(gblocks(i),:),'-o','Color',col(i,:),'LineWidth',1,'MarkerFaceColor',col(i,:))
        for i = 1:Nblock
            plot(f_x,delay.x_xAll(:,i,j),'-o','Color',col(i,:),'LineWidth',1,'MarkerFaceColor',col(i,:))
        end
        set(gca,'box','off','TickDir','out')
        title('X frequencies')
        if j == 1
            ylabel('Response Lag (ms)')
        end
        axis([0 2 0 900])
        
        subplot(2,Nsubj,j+Nsubj); hold on
%         plot(f_y,delay.y_y(gblocks(i),:),'-o','Color',col(i,:),'LineWidth',1,'MarkerFaceColor',col(i,:))
        for i = 1:Nblock
            plot(f_y,delay.y_yAll(:,i,j),'-o','Color',col(i,:),'LineWidth',1,'MarkerFaceColor',col(i,:))
        end
        set(gca,'box','off','TickDir','out')
        title('Y frequencies')
        if j == 1
            ylabel('Response Lag (ms)')
        end
        xlabel('Frequency (Hz)')
        axis([0 2 0 900])
    end
    legend(graph_name{gblocks},'Location','northeast')
end