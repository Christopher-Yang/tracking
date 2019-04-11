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
    names = {' (Rotation)',' (Mirror Reversal)'};
    
    H = gobjects(length(gblocks)+1,1);
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

        names1 = {'x','y','dx','dy'};
        names2 = {'bootx','booty','bootdx','bootdy'};
        names3 = {'errx','erry','errdx','errdy'};
        names4 = {'x_all','y_all','dx_all','dy_all'};
            
        for i = 1:Nsubj
            for j = 1:length(gblocks)
                a.x_all(:,j,i) = d.(subj{i}).(block_name{gblocks(j)}).cursor.x_fft.amplitude; % puts each subject's amplitude spectrums into data structure
                a.y_all(:,j,i) = d.(subj{i}).(block_name{gblocks(j)}).cursor.y_fft.amplitude;
                r.x(:,j,i) = d.(subj{i}).(block_name{gblocks(j)}).cursor.x_fft.fft(1:length(xAxis)); % for averaging together each subject's ffts
                r.y(:,j,i) = d.(subj{i}).(block_name{gblocks(j)}).cursor.y_fft.fft(1:length(xAxis));
            end
        end
        
        for i = 1:2
            r.(names1{i}) = mean(r.(names1{i}),3);
            a.(names1{i}) = abs(r.(names1{i})/length(r.(names1{i})));
            a.(names1{i})(2:end-1,:,:) = 2*a.(names1{i})(2:end-1,:,:);
        end
        
        a.x_all = permute(a.x_all,[1 3 2]);
        a.y_all = permute(a.y_all,[1 3 2]);
        % could alternatively not use abs() below to get a more accurate
        % measurement of aftereffects
        a.dx_all = abs(a.x_all(sort([ix;iy]),:,4)./repmat(sort([ampX ampY],'descend')',[1 Nsubj]) - a.x_all(sort([ix;iy]),:,1)./repmat(sort([ampX ampY],'descend')',[1 Nsubj]));
        a.dy_all = abs(a.y_all(sort([ix;iy]),:,4)./repmat(sort([ampX ampY],'descend')',[1 Nsubj]) - a.y_all(sort([ix;iy]),:,1)./repmat(sort([ampX ampY],'descend')',[1 Nsubj]));
        a.dx = mean(a.dx_all,2);
        a.dy = mean(a.dy_all,2);
%         diffX = a.x_all(sort([ix;iy]),:,4)./repmat(sort([ampX ampY],'descend')',[1 Nsubj]) - a.x_all(sort([ix;iy]),:,1)./repmat(sort([ampX ampY],'descend')',[1 Nsubj]);
%         diffY = a.y_all(sort([ix;iy]),:,4)./repmat(sort([ampX ampY],'descend')',[1 Nsubj]) - a.y_all(sort([ix;iy]),:,1)./repmat(sort([ampX ampY],'descend')',[1 Nsubj]);
%         a.dx_on_all = diffX(1:2:end,:);
%         a.dx_off_all = diffX(2:2:end,:);
%         a.dy_on_all = diffY(2:2:end,:);
%         a.dy_off_all = diffY(1:2:end,:);
%         diffX = mean(diffX,2);
%         diffY = mean(diffY,2);
%         a.dx_on = diffX(1:2:end);
%         a.dx_off = diffX(2:2:end);
%         a.dy_on = diffY(2:2:end);
%         a.dy_off = diffY(1:2:end);
        
        x(:,:,k) = a.dx_all';
        y(:,:,k) = a.dy_all';
        
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
            if i == 1 || i == length(gblocks)
                % plot mean with error bars; this doesn't work if averaging
                % is done in complex domain, only if individual amplitude
                % spectrums are averaged together
%                 s = shadedErrorBar(xAxis, a.x(:,i),a.errx(:,:,i));
%                 editErrorBar(s,col(3,:),0.25);
                plot(xAxis,a.x(:,i),'Color',col(3,:),'LineWidth',lw)
                plot(xAxis,a.x_all(:,:,i),'Color',[col(3,:) 0.15],'LineWidth',0.25)
                if i == length(gblocks)
                    plot(xAxis,a.x(:,1),'--','Color',col(5,:),'LineWidth',lw)
                end
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
            title([graph_name{gblocks(i)}, names{k}],'FontSize',titleSize)
            if k == 2
                legend({'X target','Y target','Hand','Baseline'})
                legend boxoff
            end
%             pbaspect([1 1 1])
            
            subplot(2,2,k+2)
            plot(freqsX, ampX,'o','LineWidth',lw,'Color',col(1,:),'MarkerFaceColor',col(1,:),'MarkerEdgeColor','none')
            hold on
            plot(freqsY, ampY,'o','LineWidth',lw,'Color',col(2,:),'MarkerFaceColor',col(2,:),'MarkerEdgeColor','none')
            if i == 1 || i == length(gblocks)
%                 s = shadedErrorBar(xAxis, a.y(:,i),a.erry(:,:,i),'lineProps','-b');
%                 editErrorBar(s,col(3,:),0.25);
                plot(xAxis,a.y(:,i),'LineWidth',lw,'Color',col(3,:),'MarkerFaceColor',col(3,:))
                plot(xAxis,a.y_all(:,:,i),'Color',[col(3,:) 0.15],'LineWidth',0.25)
                if i == length(gblocks)
                    plot(xAxis,a.y(:,1),'--','Color',col(5,:),'LineWidth',lw)
                end
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
        end
        
        figure(H(end))
        subplot(1,2,k)
        s = shadedErrorBar(sort([freqsX freqsY]),100*a.dx,100*a.errdx);
        editErrorBar(s,col(3,:),0.25);
        hold on
        s = shadedErrorBar(sort([freqsX freqsY]),100*a.dy,100*a.errdy);
        editErrorBar(s,col(4,:),0.25);
        if k == 1
            title('Rotation')
        else
            title('Mirror Reversal')
        end
        set(gca,'Xscale','log','TickDir','out')
        axis([0.1 2.3 0 60])
        xlabel('Frequency(Hz)')
        ylabel('|Post - Baseline| (%)')
        pbaspect([1 1 1])
    end
    
% %     comp.x_on = [x_on(:,:,1); x_on(:,:,2)];
% %     comp.x_off = [x_off(:,:,1); x_off(:,:,2)];
%     comp = [x_on(:,:,1); x_on(:,:,2)];
%     comp = cat(3, comp, [x_off(:,:,1); x_off(:,:,2)]);
% %     comp.y_on = [y_on(:,:,1); y_on(:,:,2)];
% %     comp.y_off = [y_off(:,:,1); y_off(:,:,2)];
%     comp2 = [repelem(0,10)'; repelem(1,10)'];
%     within = {'frequency','axis'};
%     between = {'perturbation'};
%     [tbl,rm] = simple_mixed_anova(comp, comp2, within, between)
    
end