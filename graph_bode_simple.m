function graph_bode_simple(data, graph_name, gblocks, output)
    
    col = lines;
    col = col(1:7,:);
    dat = data{end}.(output);
    
    names = {'x_x','x_y','y_y','y_x'};
    for i = 1:length(names)
        err.amp.(names{i}) = permute(dat.(names{i}).amp_err,[3 2 1]);
        err.phase.(names{i}) = permute(dat.(names{i}).phase_err,[3 2 1]);
    end
    
    f_x = dat.x_x.freqs;
    f_y = dat.y_y.freqs; 
    
    gains = [0.01 0.1 0.25 0.5 1];
    Yt = 20*log10(gains);
    Ytlab = num2cell(gains);
    legend_size = 12;
    axis_size = 12;
    label_size = 12;
    title_size = 20;
    line_size = 1.25;
    titles = {'Rotation', 'Mirror Reversal'};
    
%% graph x_x and y_y Bode plots
    figure
    subplot(2,2,1); hold on%x input, x output
    for j = 1:length(gblocks)
        s = shadedErrorBar(f_x,20*log10(dat.x_x.amplitude(gblocks(j),:)),err.amp.x_x(:,:,gblocks(j)),'lineProps','-o');
        editErrorBar(s,col(j,:),1);
    end
    set(gca,'Xscale', 'log', 'Xgrid', 'on', 'Ytick', Yt,'Yticklabel',Ytlab,'TickDir','out','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','LineWidth',1)
    title(['X ',output,' movements'])
    ylabel('Gain (cm/cm)')
    xlim([0.09 2.3])

    subplot(2,2,3); hold on
    for j = 1:length(gblocks)
        s = shadedErrorBar(f_x,dat.x_x.phase(gblocks(j),:)*(180/pi),err.phase.x_x(:,:,gblocks(j)),'lineProps','-o');
        editErrorBar(s,col(j,:),1);
    end
    set(gca, 'Xscale', 'log','Xgrid','on','TickDir','out','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','LineWidth',1)
    xlabel('Frequency (Hz)')
    ylabel(['Phase (',char(176),')'])
    xlim([0.09 2.3])

    subplot(2,2,2); hold on %y input, y output
    for j = 1:length(gblocks)
        s = shadedErrorBar(f_y,20*log10(dat.y_y.amplitude(gblocks(j),:)),err.amp.y_y(:,:,gblocks(j)),'lineProps','-o');
        editErrorBar(s,col(j,:),1);
    end
    set(gca,'Xscale','log','Xgrid','on','Ytick',Yt,'Yticklabel',Ytlab,'TickDir','out','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','LineWidth',1)
    title(['Y ',output,' movements'])
    xlim([0.09 2.3])

    subplot(2,2,4); hold on
    for j = 1:length(gblocks)
        s = shadedErrorBar(f_y,dat.y_y.phase(gblocks(j),:)*(180/pi),err.phase.y_y(:,:,gblocks(j)),'lineProps','-o');
        editErrorBar(s,col(j,:),1);
    end
    set(gca, 'Xscale', 'log','Xgrid','on','TickDir','out','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','LineWidth',1)
    xlabel('Frequency (Hz)')
    xlim([0.09 2.3])
    legend(graph_name{gblocks})

%% graph Bode plots with X and Y flipped in the graphs
    figure
    subplot(2,2,1); hold on % x input, x output
    for j = 1:length(gblocks)
        if gblocks(j) == 1
            s = shadedErrorBar(f_x,20*log10(dat.x_x.amplitude(gblocks(j),:)),err.amp.x_x(:,:,gblocks(j)),'lineProps','-o');
            editErrorBar(s,col(j,:),1);
        else
            s = shadedErrorBar(f_y,20*log10(dat.y_x.amplitude(gblocks(j),:)),err.amp.y_x(:,:,gblocks(j)),'lineProps','-o');
            editErrorBar(s,col(j,:),1);
        end
    end
    set(gca,'Xscale', 'log', 'Xgrid', 'on', 'Ytick', Yt,'Yticklabel',Ytlab,'TickDir','out','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','LineWidth',1);
    title(['X ',output,' movements']);
    ylabel('Gain (cm/cm)');
    xlim([0.09 2.3]);

    subplot(2,2,3); hold on
    for j = 1:length(gblocks)
        if gblocks(j) == 1
            s = shadedErrorBar(f_x,dat.x_x.phase(gblocks(j),:)*(180/pi),err.phase.x_x(:,:,gblocks(j)),'lineProps','-o');
            editErrorBar(s,col(j,:),1);
        else
            s = shadedErrorBar(f_y,dat.y_x.phase(gblocks(j),:)*(180/pi),err.phase.y_x(:,:,gblocks(j)),'lineProps','-o');
            editErrorBar(s,col(j,:),1);
        end
    end
    set(gca, 'Xscale', 'log','Xgrid','on','TickDir','out','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','LineWidth',1)
    xlabel('Frequency (Hz)','Fontsize',label_size)
    ylabel(['Phase (',char(176),')'],'Fontsize',label_size)
    xlim([0.09 2.3])

    subplot(2,2,2); hold on %y input, y output
    for j = 1:length(gblocks)
        if gblocks(j) == 1 || gblocks(j) == 6
            s = shadedErrorBar(f_y,20*log10(dat.y_y.amplitude(gblocks(j),:)),err.amp.y_y(:,:,gblocks(j)),'lineProps','-o');
            editErrorBar(s,col(j,:),1);
        else
            s = shadedErrorBar(f_x,20*log10(dat.x_y.amplitude(gblocks(j),:)),err.amp.x_y(:,:,gblocks(j)),'lineProps','-o');
            editErrorBar(s,col(j,:),1);
        end
    end
    set(gca,'Xscale', 'log', 'Fontsize', axis_size, 'Xgrid', 'on', 'Ytick', Yt,'Yticklabel',Ytlab,'TickDir','out','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','LineWidth',1);
    title(['Y ',output,' movements'])
    xlim([0.09 2.3]);

    subplot(2,2,4); hold on
    for j = 1:length(gblocks)
        if gblocks(j) == 1
            s = shadedErrorBar(f_y,dat.y_y.phase(gblocks(j),:)*(180/pi),err.phase.y_y(:,:,gblocks(j)),'lineProps','-o');
            editErrorBar(s,col(j,:),1);
        else
            s = shadedErrorBar(f_x,dat.x_y.phase(gblocks(j),:)*(180/pi),err.phase.x_y(:,:,gblocks(j)),'lineProps','-o');
            editErrorBar(s,col(j,:),1);
        end
    end
    set(gca, 'Xscale', 'log','Xgrid','on','TickDir','out','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','LineWidth',1);
    xlabel('Frequency (Hz)');
    xlim([0.09 2.3]);
    legend(graph_name{gblocks});
end