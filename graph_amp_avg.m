function graph_amp_avg(data,groups,block_name,gblocks,graph_name)
    
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
    
    for k = 1:length(groups)
        d = data.(groups{k});
        subj = fieldnames(d);
        subj(end) = [];
        Nsubj = length(subj);
        freqsX = d.avg.x_x.freqs;
        freqsY = d.avg.y_y.freqs;
        ampX = d.avg.x_x.amp;
        ampY = d.avg.y_y.amp;
        xAxis = d.avg.x_x.x_axis;
        ix = d.avg.x_x.index;
        iy = d.avg.y_y.index;
        
        g1 = NaN(length(xAxis),length(gblocks),Nsubj);
        g2 = NaN(length(xAxis),1000,length(gblocks));
        g3 = NaN(length(freqsX),1000);
        amps.(groups{k}).x = g1;
        amps.(groups{k}).y = g1;
        a = amps.(groups{k});
        r = amps.(groups{k});
        a.bootx = g2;
        a.booty = g2;
        a.bootdx_on = g3;
        a.bootdx_off = g3;
        a.bootdy_on = g3;
        a.bootdy_off = g3;

        names1 = {'x','y'};
        names2 = {'bootx','booty'};
        names3 = {'errx','erry'};
        names4 = {'x_all','y_all'};
            
        for i = 1:Nsubj
            for j = 1:length(gblocks)
                a.x_all(:,j,i) = d.(subj{i}).(block_name{gblocks(j)}).cursor.x_fft.amplitude; % puts each subject's amplitude spectrums into data structure
                a.y_all(:,j,i) = d.(subj{i}).(block_name{gblocks(j)}).cursor.y_fft.amplitude;
            end
        end
        
        a.x = mean(a.x_all,3);
        a.y = mean(a.y_all,3);
        
        a.x_all = permute(a.x_all,[1 3 2]);
        a.y_all = permute(a.y_all,[1 3 2]);
        
        % bootstrap data
        for j = 1:length(names1)
            for i = 1:1000
                a.(names2{j})(:,i,:) = mean(datasample(a.(names4{j}),Nsubj,2),2);
            end
            a.(names2{j}) = sort(a.(names2{j}),2);
            a.(names3{j})(:,:,1) = a.(names1{j}) - squeeze(a.(names2{j})(:,26,:)); 
            a.(names3{j})(:,:,2) = squeeze(a.(names2{j})(:,975,:)) - a.(names1{j});
            a.(names3{j}) = permute(a.(names3{j}),[3 1 2]);
        end
        
        for i = 1:length(gblocks)
            figure(H(i))
            subplot(2,2,k)
            plot(freqsX,ampX,'o','Color',col(1,:),'LineWidth',lw,'MarkerFaceColor',col(1,:),'MarkerEdgeColor','none')
            hold on
            plot(freqsY,ampY,'o','Color',col(2,:),'LineWidth',lw,'MarkerFaceColor',col(2,:),'MarkerEdgeColor','none')
            plot([50 51], [50 51],'Color',col(3,:),'LineWidth',lw)
            plot([50 51], [50 51],'--','Color',col(5,:),'LineWidth',lw)
                % plot mean with error bars; doesn't work if amplitude
                % spectra are averaged in the complex domain
%                 s = shadedErrorBar(xAxis, a.x(:,i),a.errx(:,:,i));
%                 editErrorBar(s,col(3,:),0.25);
            plot(xAxis,a.x(:,i),'Color',col(3,:),'LineWidth',lw)
            plot(xAxis,a.x_all(:,:,i),'Color',[col(3,:) 0.15],'LineWidth',0.25)
            if i ~= 1
                plot(xAxis,a.x(:,1),'--','Color',col(5,:),'LineWidth',lw)
            end
            plot(freqsX,a.x(ix,i),'o','Color',col(1,:),'LineWidth',lw,'MarkerFaceColor',col(1,:),'MarkerEdgeColor','none')

            set(gca,'Xtick',0:1:2,'Ytick',0:0.01:0.03,'box','off','LineWidth',1,'TickDir','out','FontSize',fontSize,'XMinorTick','on')
            ax = gca;
            ax.XAxis.MinorTickValues = 0:0.25:2.25;
            axis([0 2.3 0 ampX(1)+.01])
            ylabel('X Amplitude (m)')
            title(graph_name{gblocks(i)},'FontSize',titleSize)
            if k == 2
                legend({'X target','Y target','Hand','Baseline'})
                legend boxoff
            end
%             pbaspect([1 1 1])
            
            subplot(2,2,k+2)
            plot(freqsX, ampX,'o','LineWidth',lw,'Color',col(1,:),'MarkerFaceColor',col(1,:),'MarkerEdgeColor','none')
            hold on
            plot(freqsY, ampY,'o','LineWidth',lw,'Color',col(2,:),'MarkerFaceColor',col(2,:),'MarkerEdgeColor','none')

%                 s = shadedErrorBar(xAxis, a.y(:,i),a.erry(:,:,i),'lineProps','-b');
%                 editErrorBar(s,col(3,:),0.25);
            plot(xAxis,a.y(:,i),'LineWidth',lw,'Color',col(3,:),'MarkerFaceColor',col(3,:))
            plot(xAxis,a.y_all(:,:,i),'Color',[col(3,:) 0.15],'LineWidth',0.25)
            if i ~= 1
                plot(xAxis,a.y(:,1),'--','Color',col(5,:),'LineWidth',lw)
            end
            plot(freqsY,a.y(iy,i),'o','LineWidth',1.5,'Color',col(2,:),'MarkerFaceColor',col(2,:),'MarkerEdgeColor','none')

            set(gca,'Xtick',0:1:2,'Ytick',0:0.01:0.03,'box','off','LineWidth',1,'TickDir','out','FontSize',fontSize,'XMinorTick','on')
            ax = gca;
            ax.XAxis.MinorTickValues = 0:0.25:2.25;
            axis([0 2.3 0 ampY(1)+.01])
            ylabel('Y Amplitude (m)')
            xlabel('Frequency (Hz)')
%             pbaspect([1 1 1])
        end
    end
end