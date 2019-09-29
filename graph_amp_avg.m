function graph_amp_avg(data,block_name,gblocks,graph_name,output)
    
    col = [30 144 255
           0 0 139
           255 182 193]./255;
    lw = 1;
    
    H = gobjects(length(gblocks),1);
    for i = 1:length(H)
        H(i) = figure;
    end
    
    Nsubj = length(data)-1;
    freqsX = data{end}.freqX;
    freqsY = data{end}.freqY;
    ampX = data{end}.ampX;
    ampY = data{end}.ampY;
    xAxis = data{end}.x_axis;
    ix = data{end}.(output).x_x.index;
    iy = data{end}.(output).y_y.index;
    Nsteps = length(data{1}.(block_name{1}).(output).x_pos);
    
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
            dat = data{i}.(block_name{gblocks(j)});
            pX = dat.(output).x_fft;
            pY = dat.(output).y_fft;
            a.x_all(:,j,i) = squeeze(abs(mean(pX,2)));
            a.y_all(:,j,i) = squeeze(abs(mean(pY,2)));
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
        subplot(2,1,1); hold on
        plot(freqsX,ampX,'ko','LineWidth',lw,'MarkerFaceColor',[0 0 0])
        plot(freqsY,ampY,'ko','LineWidth',lw)
        % plot mean with error bars; doesn't work if amplitude spectra are averaged in the complex domain
%         s = shadedErrorBar(xAxis, a.x(:,i),a.errx(:,:,i));
%         editErrorBar(s,col(3,:),0.25);
        plot(xAxis,a.x(:,i),'Color',col(1,:),'LineWidth',lw)
        plot(freqsX,a.x(ix,i),'ko-','LineWidth',lw,'MarkerFaceColor',[0 0 0])
        plot(freqsY,a.x(iy,i),'ko-','LineWidth',lw)
        set(gca,'Xtick',0:1:2,'Ytick',0:0.01:0.03,'box','off','LineWidth',1,'TickDir','out','XMinorTick','on')
        ax = gca;
        ax.XAxis.MinorTickValues = 0:0.25:2.25;
        axis([0 2.3 0 0.025])
        ylabel('X Amplitude (m)')
        title(graph_name{gblocks(i)})
        legend({'X target','Y target','Hand'})
        legend boxoff
        
        subplot(2,1,2); hold on
        plot(freqsX, ampX,'ko','LineWidth',lw,'MarkerFaceColor',[0 0 0])
        plot(freqsY, ampY,'ko','LineWidth',lw)
%         s = shadedErrorBar(xAxis, a.y(:,i),a.erry(:,:,i),'lineProps','-b');
%         editErrorBar(s,col(3,:),0.25);
        plot(xAxis,a.y(:,i),'LineWidth',lw,'Color',col(1,:))
        plot(freqsY,a.y(iy,i),'ko-','LineWidth',lw)
        plot(freqsX,a.y(ix,i),'ko-','LineWidth',lw,'MarkerFaceColor',[0 0 0])
        set(gca,'Xtick',0:1:2,'Ytick',0:0.01:0.03,'box','off','LineWidth',1,'TickDir','out','XMinorTick','on')
        ax = gca;
        ax.XAxis.MinorTickValues = 0:0.25:2.25;
        axis([0 2.3 0 0.025])
        ylabel('Y Amplitude (m)')
        xlabel('Frequency (Hz)')
    end
end