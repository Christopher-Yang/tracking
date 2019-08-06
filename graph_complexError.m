function graph_complexError(data, graph_name, gblocks)
    
    group_names = {'Rotation','Mirror Reversal'};
    f_x = data.(groups{1}).avg.x_x.freqs;
    f_y = data.(groups{1}).avg.y_y.freqs;
    a = figure('Name','Complex Tracking Error','NumberTitle','off');
    b = figure('Name','Difference Between Ideal and Observed Tracking Error','NumberTitle','off');
    line_width = 1.5;
    font_size = 12;
    
    for i = 1:length(groups)
        x_x = NaN(size(data.(groups{i}).avg.x_x.fft));
        y_y = NaN(size(data.(groups{i}).avg.y_y.fft));
        x_x_opt = NaN(size(data.(groups{i}).avg.x_x.optimal_fft));
        y_y_opt = NaN(size(data.(groups{i}).avg.y_y.optimal_fft));
        for j = 1:numel(x_x)
            x_x(j) = norm(1 - data.(groups{i}).avg.x_x.fft(j));
            y_y(j) = norm(1 - data.(groups{i}).avg.y_y.fft(j));
            x_x_opt(j) = norm(1 - data.(groups{i}).avg.x_x.optimal_fft(j));
            y_y_opt(j) = norm(1 - data.(groups{i}).avg.y_y.optimal_fft(j));
        end
        figure(a)
        subplot(2,2,i)
%         semilogx(f_x,x_x(gblocks,:)','-o','LineWidth',line_width)
        semilogx(f_x,x_x(1,:)','-o','LineWidth',line_width)
        hold on
        semilogx(f_y,y_y([2 5],:)','-o','LineWidth',line_width)
        semilogx(f_x,x_x(6,:)','-o','LineWidth',line_width)
        set(gca,'ColorOrderIndex',1,'box','off','xgrid','on','LineWidth',1,'FontSize',font_size)
%         semilogx(f_x,x_x_opt','--o')
        title(group_names{i},'FontSize',font_size)
        xlabel('Frequency (Hz)')
        ylabel('X Tracking Error')
        axis([0.09 2.3 0 1.8])
        
        subplot(2,2,i+2)
%         semilogx(f_y,y_y(gblocks,:)','-o','LineWidth',line_width)
        semilogx(f_y,y_y(1,:)','-o','LineWidth',line_width)
        hold on
        semilogx(f_x,x_x([2 5],:)','-o','LineWidth',line_width)
        semilogx(f_y,y_y(6,:)','-o','LineWidth',line_width)
        set(gca,'ColorOrderIndex',1,'box','off','xgrid','on','LineWidth',1,'FontSize',font_size);
%         semilogx(f_y,y_y_opt','--o');
        xlabel('Frequency (Hz)'); 
        ylabel('Y Tracking Error');
        axis([0.09 2.3 0 1.8]);
        legend(graph_name,'Position',[0.87 0.4 0.1 0.2],'FontSize',font_size);
        
        figure(b)
        subplot(2,2,i)
%         semilogx(f_x,x_x(gblocks,:)' - x_x_opt(gblocks,:)','-o','LineWidth',3)
        semilogx(f_x,x_x(1,:)' - x_x_opt(1,:)','-o','LineWidth',line_width)
        hold on
        semilogx(f_y,y_y([2 5],:)' - y_y_opt([2 5],:)','-o','LineWidth',line_width)
        semilogx(f_x,x_x(6,:)' - x_x_opt(6,:)','-o','LineWidth',line_width)
        set(gca,'box','off','xgrid','on','LineWidth',1,'FontSize',font_size)
        title(group_names{i},'FontSize',font_size)
        xlabel('Frequency (Hz)')
        ylabel('||X data - X optimal||','FontSize',font_size)
        axis([0.09 2.3 0 1.8])
        
        subplot(2,2,i+2)
%         semilogx(f_y,y_y(gblocks,:)' - y_y_opt(gblocks,:)','-o','LineWidth',3)
        semilogx(f_y,y_y(1,:)' - y_y_opt(1,:)','-o','LineWidth',line_width)
        hold on
        semilogx(f_x,x_x([2 5],:)' - x_x_opt([2 5],:)','-o','LineWidth',line_width)
        semilogx(f_y,y_y(6,:)' - y_y_opt(6,:)','-o','LineWidth',line_width)
        set(gca,'box','off','xgrid','on','LineWidth',1,'FontSize',font_size)
        xlabel('Frequency (Hz)')
        ylabel('||Y data - Y optimal||','FontSize',font_size)
        axis([0.09 2.3 0 1.8])
        legend(graph_name,'Position',[0.87 0.4 0.1 0.2],'FontSize',font_size)
    end
end