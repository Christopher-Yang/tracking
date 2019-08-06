function graph_amp_avg(data,block_name,gblocks,graph_name)
    
    col = [255 69 0
           153 50 204
           0 0 139
           255 182 193
           199 178 153
           128 0 0]./255;
    fontSize = 12;
    titleSize = 14;
    lw = 1;
    
    H = gobjects(length(gblocks),1);
    for i = 1:length(H)
        H(i) = figure;
    end
    
%     d = data;
%     subj = fieldnames(d);
%     subj(end) = [];
    Nsubj = length(data)-1;
    freqsX = data{end}.x_x.freqs;
    freqsY = data{end}.y_y.freqs;
    ampX = data{end}.x_x.amp;
    ampY = data{end}.y_y.amp;
    xAxis = data{end}.x_x.x_axis;
    ix = data{end}.x_x.index;
    iy = data{end}.y_y.index;
    Nsteps = length(data{1}.(block_name{1}).cursor.x_pos);
    
    g1 = NaN(length(xAxis),length(gblocks),Nsubj);
    g2 = NaN(length(xAxis),1000,length(gblocks));
    amps.x = g1;
    amps.y = g1;
    a = amps;
    
    names1 = {'x','y'};
    names2 = {'bootx','booty'};
    names3 = {'errx','erry'};
    names4 = {'x_all','y_all'};
    
    for i = 1:Nsubj
        for j = 1:length(gblocks)
            a.x_all(:,j,i) = data{i}.(block_name{gblocks(j)}).cursor.x_fft.amplitude; % puts each subject's amplitude spectrums into data structure
            a.y_all(:,j,i) = data{i}.(block_name{gblocks(j)}).cursor.y_fft.amplitude;
        end
    end
    
    a.x = mean(a.x_all,3);
    a.y = mean(a.y_all,3);
    
    a.x_all = permute(a.x_all,[1 3 2]);
    a.y_all = permute(a.y_all,[1 3 2]);
    
    % bootstrap data
%     for j = 1:length(names1)
%         for i = 1:1000
%             a.(names2{j})(:,i,:) = mean(datasample(a.(names4{j}),Nsubj,2),2);
%         end
%         a.(names2{j}) = sort(a.(names2{j}),2);
%         a.(names3{j})(:,:,1) = a.(names1{j}) - squeeze(a.(names2{j})(:,26,:));
%         a.(names3{j})(:,:,2) = squeeze(a.(names2{j})(:,975,:)) - a.(names1{j});
%         a.(names3{j}) = permute(a.(names3{j}),[3 1 2]);
%     end
    
    for i = 1:length(gblocks)
        figure(H(i))
        subplot(2,1,1)
        plot(freqsX,ampX,'o','Color',col(1,:),'LineWidth',lw,'MarkerFaceColor',col(1,:),'MarkerEdgeColor','none')
        hold on
        plot(freqsY,ampY,'o','Color',col(2,:),'LineWidth',lw,'MarkerFaceColor',col(2,:),'MarkerEdgeColor','none')
        plot([50 51], [50 51],'Color',col(3,:),'LineWidth',lw)
        plot([50 51], [50 51],'--','Color',col(5,:),'LineWidth',lw)
        if i == 1
            % plot mean with error bars; doesn't work if amplitude
            % spectra are averaged in the complex domain
            %                 s = shadedErrorBar(xAxis, a.x(:,i),a.errx(:,:,i));
            %                 editErrorBar(s,col(3,:),0.25);
            plot(xAxis,a.x(:,i),'Color',col(3,:),'LineWidth',lw)
            plot(xAxis,a.x_all(:,:,i),'Color',[col(3,:) 0.15],'LineWidth',0.25)
            plot(freqsX,a.x(ix,i),'o','Color',col(1,:),'LineWidth',lw,'MarkerFaceColor',col(1,:),'MarkerEdgeColor','none')
        else
            %                 s = shadedErrorBar(xAxis, a.y(:,i),a.erry(:,:,i));
            %                 editErrorBar(s,col(3,:),0.25);
            plot(xAxis,a.y(:,i),'Color',col(3,:),'LineWidth',lw)
            plot(xAxis,a.y_all(:,:,i),'Color',[col(3,:) 0.15],'LineWidth',0.25)
            plot(xAxis,a.x(:,1),'--','Color',col(5,:),'LineWidth',lw)
            plot(freqsY,a.y(iy,i),'o','Color',col(2,:),'LineWidth',lw,'MarkerFaceColor',col(2,:),'MarkerEdgeColor','none')
        end
        set(gca,'Xtick',0:1:2,'Ytick',0:0.01:0.03,'box','off','LineWidth',1,'TickDir','out','FontSize',fontSize,'XMinorTick','on')
        ax = gca;
        ax.XAxis.MinorTickValues = 0:0.25:2.25;
        axis([0 2.3 0 ampX(1)+.01])
        ylabel('X Amplitude (m)')
        title(graph_name{gblocks(i)},'FontSize',titleSize)
        legend({'X target','Y target','Hand','Baseline'})
        legend boxoff
        
        subplot(2,1,2)
        plot(freqsX, ampX,'o','LineWidth',lw,'Color',col(1,:),'MarkerFaceColor',col(1,:),'MarkerEdgeColor','none')
        hold on
        plot(freqsY, ampY,'o','LineWidth',lw,'Color',col(2,:),'MarkerFaceColor',col(2,:),'MarkerEdgeColor','none')
        if i == 1
            %                 s = shadedErrorBar(xAxis, a.y(:,i),a.erry(:,:,i),'lineProps','-b');
            %                 editErrorBar(s,col(3,:),0.25);
            plot(xAxis,a.y(:,i),'LineWidth',lw,'Color',col(3,:),'MarkerFaceColor',col(3,:))
            plot(xAxis,a.y_all(:,:,i),'Color',[col(3,:) 0.15],'LineWidth',0.25)
            plot(freqsY,a.y(iy,i),'o','LineWidth',1.5,'Color',col(2,:),'MarkerFaceColor',col(2,:),'MarkerEdgeColor','none')
        else
            %                 s = shadedErrorBar(xAxis, a.x(:,i),a.errx(:,:,i),'lineProps','-b');
            %                 editErrorBar(s,col(3,:),0.25);
            plot(xAxis,a.x(:,i),'LineWidth',lw,'Color',col(3,:),'MarkerFaceColor',col(3,:))
            plot(xAxis,a.x_all(:,:,i),'Color',[col(3,:) 0.15],'LineWidth',0.25)
            plot(xAxis,a.y(:,1),'--','Color',col(5,:),'LineWidth',lw)
            plot(freqsX,a.x(ix,i),'o','LineWidth',lw,'Color',col(1,:),'MarkerFaceColor',col(1,:),'MarkerEdgeColor','none')
        end
        set(gca,'Xtick',0:1:2,'Ytick',0:0.01:0.03,'box','off','LineWidth',1,'TickDir','out','FontSize',fontSize,'XMinorTick','on')
        ax = gca;
        ax.XAxis.MinorTickValues = 0:0.25:2.25;
        axis([0 2.3 0 ampY(1)+.01])
        ylabel('Y Amplitude (m)')
        xlabel('Frequency (Hz)')
        %             pbaspect([1 1 1])
    end
end