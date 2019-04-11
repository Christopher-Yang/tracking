function graph_bode_simple(data, graph_name, groups, gblocks)
    
%     err = struct('rot',NaN(6,7,2),'mir',NaN(6,7,2));
%     for i = 1:2
%         err.(groups{i}).x_x(:,:,1) = data.(groups{i}).avg.x_x.d.amp_err(:,:,1) - 20*log10(data.(groups{i}).avg.x_x.d.amplitude);
%         err.(groups{i}).x_x(:,:,2) = 20*log10(data.(groups{i}).avg.x_x.d.amplitude) - data.(groups{i}).avg.x_x.d.amp_err(:,:,2);
%         err.(groups{i}).y_y(:,:,1) = data.(groups{i}).avg.y_y.d.amp_err(:,:,1) - 20*log10(data.(groups{i}).avg.y_y.d.amplitude);
%         err.(groups{i}).y_y(:,:,2) = 20*log10(data.(groups{i}).avg.y_y.d.amplitude) - data.(groups{i}).avg.y_y.d.amp_err(:,:,2);
%     end
    
    col = lines;
    col = col(1:7,:);
    
    for i = 1:length(groups)
        err.amp.(groups{i}).x_x = permute(data.(groups{i}).avg.x_x.amp_err,[3 2 1]);
        err.amp.(groups{i}).y_y = permute(data.(groups{i}).avg.y_y.amp_err,[3 2 1]);
        err.phase.(groups{i}).x_x = permute(data.(groups{i}).avg.x_x.phase_err,[3 2 1]);
        err.phase.(groups{i}).y_y = permute(data.(groups{i}).avg.y_y.phase_err,[3 2 1]);
    end
    
    f_x = data.(groups{1}).avg.x_x.freqs;
    f_y = data.(groups{1}).avg.y_y.freqs; 
    
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
    m = figure('Name','Bode (X)','NumberTitle','off','DefaultTextFontName','Arial','DefaultAxesFontName','Arial');
    n = figure('Name','Bode (Y)','NumberTitle','off','DefaultTextFontName','Arial','DefaultAxesFontName','Arial');
    for i = 1:length(groups)
        figure(m);
        subplot(2,2,i); %x input, x output
        for j = 1:length(gblocks)
            s = shadedErrorBar(f_x,20*log10(data.(groups{i}).avg.x_x.amplitude(gblocks(j),:)),err.amp.(groups{i}).x_x(:,:,gblocks(j)),'lineProps','-o');
            editErrorBar(s,col(j,:),1);
            hold on
        end
        set(gca,'Xscale', 'log', 'Fontsize', axis_size, 'Xgrid', 'on', 'Ytick', Yt,'Yticklabel',Ytlab,'TickDir','out','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','LineWidth',1);
        title(titles{i},'Fontsize',title_size); 
        if i == 1
            ylabel('Gain (cm/cm)','Fontsize',label_size); 
        end
        xlim([0.09 2.3]);
        ylim([-40 2]);
        
        subplot(2,2,i+2);
        for j = 1:length(gblocks)
            s = shadedErrorBar(f_x,data.(groups{i}).avg.x_x.phase(gblocks(j),:)*(180/pi),err.phase.(groups{i}).x_x(:,:,gblocks(j)),'lineProps','-o');
            editErrorBar(s,col(j,:),1);
            hold on
        end
        set(gca, 'Xscale', 'log', 'Fontsize', axis_size,'Xgrid','on', 'Ytick',-800:100:0,'TickDir','out','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','LineWidth',1); 
        xlabel('Frequency (Hz)','Fontsize',label_size);
        if i == 1
            ylabel(['Phase (',char(176),')'],'Fontsize',label_size); 
        end
        xlim([0.09 2.3]);
        ylim([-400 0]);
        legend(graph_name{gblocks},'Position',[0.15 0.16 0.1 0.1],'FontSize',legend_size);
        
        figure(n)
        subplot(2,2,i); %y input, y output
        for j = 1:length(gblocks)
            s = shadedErrorBar(f_y,20*log10(data.(groups{i}).avg.y_y.amplitude(gblocks(j),:)),err.amp.(groups{i}).y_y(:,:,gblocks(j)),'lineProps','-o');
            editErrorBar(s,col(j,:),1);
            hold on
        end
        set(gca,'Xscale', 'log', 'Fontsize', axis_size, 'Xgrid', 'on', 'Ytick', Yt,'Yticklabel',Ytlab,'TickDir','out','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','LineWidth',1);
        title(titles{i},'Fontsize',title_size); 
        if i == 1
            ylabel('Gain (cm/cm)','Fontsize',label_size); 
        end
        xlim([0.09 2.3]);
        ylim([-40 2]);
        
        subplot(2,2,i+2);
        for j = 1:length(gblocks)
            s = shadedErrorBar(f_y,data.(groups{i}).avg.y_y.phase(gblocks(j),:)*(180/pi),err.phase.(groups{i}).y_y(:,:,gblocks(j)),'lineProps','-o');
            editErrorBar(s,col(j,:),1);
            hold on
        end
        set(gca, 'Xscale', 'log', 'Fontsize', axis_size,'Xgrid','on', 'Ytick', -800:100:0,'TickDir','out','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','LineWidth',1); 
        xlabel('Frequency (Hz)','Fontsize',label_size);
        if i == 1
            ylabel(['Phase (',char(176),')'],'Fontsize',label_size); 
        end
        xlim([0.09 2.3]);
        ylim([-400 0]);
        legend(graph_name{gblocks},'Position',[0.15 0.15 0.1 0.1],'FontSize',legend_size);
    end
%% graph hand Bode plots
%     o = figure('Name','Bode (X \rightarrow Y)','NumberTitle','off','DefaultTextFontName','Arial','DefaultAxesFontName','Arial');
%     p = figure('Name','Bode (Y \rightarrow X)','NumberTitle','off','DefaultTextFontName','Arial','DefaultAxesFontName','Arial');
%     for i = 1:length(groups)
%         figure(o); 
%         subplot(2,2,i); hold on;
%         errorbar(freqs_x(:,1), 20*log10(data.(groups{i}).avg.x_x.d.amplitude(1,:)'), data.(groups{i}).avg.x_x.d.amp_err(1,:,1)',data.(groups{i}).avg.x_x.d.amp_err(1,:,2)','LineWidth',line_size);
%         errorbar(freqs_x(:,2:3), 20*log10(data.(groups{i}).avg.x_y.d.amplitude([2 5],:)'), data.(groups{i}).avg.x_y.d.amp_err([2 5],:,1)',data.(groups{i}).avg.x_y.d.amp_err([2 5],:,2)','LineWidth',line_size);
%         errorbar(freqs_x(:,4), 20*log10(data.(groups{i}).avg.x_x.d.amplitude(6,:)'), data.(groups{i}).avg.x_x.d.amp_err(6,:,1)',data.(groups{i}).avg.x_x.d.amp_err(6,:,2)','LineWidth',line_size);
%         title(titles{i},'FontSize',title_size);
%         if i == 1
%             ylabel('Magnitude (dB)','FontSize',label_size);
%         end
%         set(gca,'Xscale', 'log', 'Fontsize', axis_size, 'Xgrid', 'on', 'Ytick', -40:10:0,'TickDir','out','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off');
%         axis([0.09 2.3 -40 1]);
%         
%         subplot(2,2,i+2); hold on;
%         errorbar(freqs_x(:,1), data.(groups{i}).avg.x_x.d.phase(1,:)'*(180/pi), data.(groups{i}).avg.x_x.d.phase_err(1,:,1)',data.(groups{i}).avg.x_x.d.phase_err(1,:,2)','LineWidth',line_size);
%         errorbar(freqs_x(:,2:3), data.(groups{i}).avg.x_y.d.phase([2 5],:)'*(180/pi), data.(groups{i}).avg.x_y.d.phase_err([2 5],:,1)',data.(groups{i}).avg.x_y.d.phase_err([2 5],:,2)','LineWidth',line_size);
%         errorbar(freqs_x(:,4), data.(groups{i}).avg.x_x.d.phase(6,:)'*(180/pi), data.(groups{i}).avg.x_x.d.phase_err(6,:,1)',data.(groups{i}).avg.x_x.d.phase_err(6,:,2)','LineWidth',line_size);
%         xlabel('Frequency (Hz)','FontSize',label_size);
%         if i == 1
%             ylabel('Phase (degrees)','FontSize',label_size); 
%         end
%         set(gca,'Xscale', 'log', 'Fontsize', axis_size, 'Xgrid', 'on', 'Ytick', -400:100:100,'TickDir','out','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off');
%         axis([0.09 2.3 -400 100]);
%         legend(graph_name,'Position',[0.15 0.61 0.1 0.1],'FontSize',legend_size);
%         
%         figure(p)
%         subplot(2,2,i); hold on;
%         errorbar(freqs_y(:,1), 20*log10(data.(groups{i}).avg.y_y.d.amplitude(1,:)'), data.(groups{i}).avg.x_x.d.amp_err(1,:,1)',data.(groups{i}).avg.x_x.d.amp_err(1,:,2)','LineWidth',line_size);
%         errorbar(freqs_y(:,2:3), 20*log10(data.(groups{i}).avg.y_x.d.amplitude([2 5],:)'), data.(groups{i}).avg.x_y.d.amp_err([2 5],:,1)',data.(groups{i}).avg.x_y.d.amp_err([2 5],:,2)','LineWidth',line_size);
%         errorbar(freqs_y(:,4), 20*log10(data.(groups{i}).avg.y_y.d.amplitude(6,:)'), data.(groups{i}).avg.x_x.d.amp_err(6,:,1)',data.(groups{i}).avg.x_x.d.amp_err(6,:,2)','LineWidth',line_size);
%         title(titles{i},'FontSize',title_size); 
%         if i == 1
%             ylabel('Magnitude (dB)','FontSize',label_size);
%         end
%         set(gca,'Xscale', 'log', 'Fontsize', axis_size, 'Xgrid', 'on', 'Ytick', -40:10:0,'TickDir','out','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off');
%         axis([0.09 2.3 -40 1]);
%         
%         subplot(2,2,i+2); hold on;
%         errorbar(freqs_y(:,1), data.(groups{i}).avg.y_y.d.phase(1,:)'*(180/pi), data.(groups{i}).avg.x_x.d.phase_err(1,:,1)',data.(groups{i}).avg.x_x.d.phase_err(1,:,2)','LineWidth',line_size);
%         errorbar(freqs_y(:,2:3), data.(groups{i}).avg.y_x.d.phase([2 5],:)'*(180/pi), data.(groups{i}).avg.x_y.d.phase_err([2 5],:,1)',data.(groups{i}).avg.x_y.d.phase_err([2 5],:,2)','LineWidth',line_size);
%         errorbar(freqs_y(:,4), data.(groups{i}).avg.y_y.d.phase(6,:)'*(180/pi), data.(groups{i}).avg.x_x.d.phase_err(6,:,1)',data.(groups{i}).avg.x_x.d.phase_err(6,:,2)','LineWidth',line_size);
%         xlabel('Frequency (Hz)','FontSize',label_size);
%         if i == 1
%             ylabel('Phase (degrees)','FontSize',label_size);
%         end
%         set(gca,'Xscale', 'log', 'Fontsize', axis_size, 'Xgrid', 'on', 'Ytick', -400:100:100,'TickDir','out','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off');
%         axis([0.09 2.3 -400 100]);
%         legend(graph_name,'Position',[0.15 0.61 0.1 0.1],'FontSize',legend_size);
%     end
    
%% graph x_y and y_x Bode plots
%     for i = 1:length(groups)
%         figure('Name',titles{i},'NumberTitle','off'); 
%         subplot(2,2,1); %x input, y output
%         errorbar(freqs_x, 20*log10(data.(groups{i}).avg.x_y.d.amplitude(gblocks,:)'), data.(groups{i}).avg.x_y.d.amp_err(gblocks,:,1)',data.(groups{i}).avg.x_y.d.amp_err(gblocks,:,2)'); 
%         title('Magnitude (X_{hand} -> Y_{cursor})'); ylabel('Magnitude (dB)'); xlabel('Frequency (Hz)');
%         grid on; set(gca, 'Xscale', 'log');
%         xlim([0.09 2.3]);
%         ylim([-60 0]);
%         
%         subplot(2,2,3);
%         errorbar(freqs_x, data.(groups{i}).avg.x_y.d.phase(gblocks,:)'*(180/pi), data.(groups{i}).avg.x_y.d.phase_err(gblocks,:,1)',data.(groups{i}).avg.x_y.d.phase_err(gblocks,:,2)');
%         title('Phase (X_{hand} -> Y_{cursor})'); ylabel('Phase (degrees)'); xlabel('Frequency (Hz)');
%         grid on; set(gca, 'Xscale', 'log');
%         xlim([0.09 2.3]);
%         ylim([-800 200]);
% 
%         subplot(2,2,2); %y input, x output
%         errorbar(freqs_y, 20*log10(data.(groups{i}).avg.y_x.d.amplitude(gblocks,:)'), data.(groups{i}).avg.y_x.d.amp_err(gblocks,:,1)',data.(groups{i}).avg.y_x.d.amp_err(gblocks,:,2)'); 
%         title('Magnitude (Y_{hand} -> X_{cursor})'); ylabel('Magnitude (dB)'); xlabel('Frequency (Hz)');
%         grid on; set(gca, 'Xscale', 'log');
%         xlim([0.09 2.3]);
%         ylim([-60 0]); 
% 
%         subplot(2,2,4);
%         errorbar(freqs_y, data.(groups{i}).avg.y_x.d.phase(gblocks,:)'*(180/pi), data.(groups{i}).avg.y_x.d.phase_err(gblocks,:,1)',data.(groups{i}).avg.y_x.d.phase_err(gblocks,:,2)');
%         title('Phase (Y_{hand} -> X_{cursor})'); ylabel('Phase (degrees)'); xlabel('Frequency (Hz)');
%         grid on; set(gca, 'Xscale', 'log');
%         xlim([0.09 2.3]);
%         ylim([-800 200]);
%         legend(graph_name,'Position',[0 0.4 0.1 0.2]);
%     end    

%% graph Bode plots with X and Y flipped in the graphs
    q = figure('Name','Bode (X flipped)','NumberTitle','off','DefaultTextFontName','Arial','DefaultAxesFontName','Arial');
    r = figure('Name','Bode (Y flipped)','NumberTitle','off','DefaultTextFontName','Arial','DefaultAxesFontName','Arial');
    for i = 1:length(groups)
        figure(q);
        subplot(2,2,i); %x input, x output
        for j = 1:length(gblocks)
            if gblocks(j) == 1 || gblocks(j) == 6
                s = shadedErrorBar(f_x,20*log10(data.(groups{i}).avg.x_x.amplitude(gblocks(j),:)),err.amp.(groups{i}).x_x(:,:,gblocks(j)),'lineProps','-o');
                editErrorBar(s,col(j,:),1);
                hold on
            else
                s = shadedErrorBar(f_y,20*log10(data.(groups{i}).avg.y_y.amplitude(gblocks(j),:)),err.amp.(groups{i}).y_y(:,:,gblocks(j)),'lineProps','-o');
                editErrorBar(s,col(j,:),1);
                hold on
            end
        end
        set(gca,'Xscale', 'log', 'Fontsize', axis_size, 'Xgrid', 'on', 'Ytick', Yt,'Yticklabel',Ytlab,'TickDir','out','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','LineWidth',1);
        title(titles{i},'Fontsize',title_size); 
        if i == 1
            ylabel('Gain (cm/cm)','Fontsize',label_size); 
        end
        xlim([0.09 2.3]);
        ylim([-40 2]);
        
        subplot(2,2,i+2);
        for j = 1:length(gblocks)
            if gblocks(j) == 1 || gblocks(j) == 6
                s = shadedErrorBar(f_x,data.(groups{i}).avg.x_x.phase(gblocks(j),:)*(180/pi),err.phase.(groups{i}).x_x(:,:,gblocks(j)),'lineProps','-o');
                editErrorBar(s,col(j,:),1);
                hold on
            else
                s = shadedErrorBar(f_y,data.(groups{i}).avg.y_y.phase(gblocks(j),:)*(180/pi),err.phase.(groups{i}).y_y(:,:,gblocks(j)),'lineProps','-o');
                editErrorBar(s,col(j,:),1);
                hold on
            end
        end
        set(gca, 'Xscale', 'log', 'Fontsize', axis_size,'Xgrid','on', 'Ytick',-800:100:0,'TickDir','out','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','LineWidth',1); 
        xlabel('Frequency (Hz)','Fontsize',label_size);
        if i == 1
            ylabel(['Phase (',char(176),')'],'Fontsize',label_size); 
        end
        xlim([0.09 2.3]);
        ylim([-400 0]);
        legend(graph_name{gblocks},'Position',[0.15 0.16 0.1 0.1],'FontSize',legend_size);
        
        figure(r)
        subplot(2,2,i); %y input, y output
        for j = 1:length(gblocks)
            if gblocks(j) == 1 || gblocks(j) == 6
                s = shadedErrorBar(f_y,20*log10(data.(groups{i}).avg.y_y.amplitude(gblocks(j),:)),err.amp.(groups{i}).y_y(:,:,gblocks(j)),'lineProps','-o');
                editErrorBar(s,col(j,:),1);
                hold on
            else
                s = shadedErrorBar(f_x,20*log10(data.(groups{i}).avg.x_x.amplitude(gblocks(j),:)),err.amp.(groups{i}).x_x(:,:,gblocks(j)),'lineProps','-o');
                editErrorBar(s,col(j,:),1);
                hold on
            end
        end
        set(gca,'Xscale', 'log', 'Fontsize', axis_size, 'Xgrid', 'on', 'Ytick', Yt,'Yticklabel',Ytlab,'TickDir','out','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','LineWidth',1);
        title(titles{i},'Fontsize',title_size); 
        if i == 1
            ylabel('Gain (cm/cm)','Fontsize',label_size); 
        end
        xlim([0.09 2.3]);
        ylim([-40 2]);
        
        subplot(2,2,i+2);
        for j = 1:length(gblocks)
            if gblocks(j) == 1 || gblocks(j) == 6
                s = shadedErrorBar(f_y,data.(groups{i}).avg.y_y.phase(gblocks(j),:)*(180/pi),err.phase.(groups{i}).y_y(:,:,gblocks(j)),'lineProps','-o');
                editErrorBar(s,col(j,:),1);
                hold on
            else
                s = shadedErrorBar(f_x,data.(groups{i}).avg.x_x.phase(gblocks(j),:)*(180/pi),err.phase.(groups{i}).x_x(:,:,gblocks(j)),'lineProps','-o');
                editErrorBar(s,col(j,:),1);
                hold on
            end
        end
        set(gca, 'Xscale', 'log', 'Fontsize', axis_size,'Xgrid','on', 'Ytick', -800:100:0,'TickDir','out','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','LineWidth',1); 
        xlabel('Frequency (Hz)','Fontsize',label_size);
        if i == 1
            ylabel(['Phase (',char(176),')'],'Fontsize',label_size); 
        end
        xlim([0.09 2.3]);
        ylim([-400 0]);
        legend(graph_name{gblocks},'Position',[0.15 0.15 0.1 0.1],'FontSize',legend_size);
    end
end