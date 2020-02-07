function graph_ampSpectra(data)
% plots the right hand's amplitude spectra

% set variables for plotting
output = 'Rhand';
groups = {'rot','mir'};
block_name = {'baseline','pert1','pert2','pert3','pert4','post'};
graph_name = {'Baseline','Early','Train2','Train3','Late','Post'};
gblocks = [1 2 5 6];
col = [30 144 255
    0 0 139
    255 182 193]./255; % colors used for plotting
lw = 1; % set width of lines
names = {' (Rotation)',' (Mirror Reversal)'}; 

for k = 1:length(groups)
    % set more variables for plotting
    d = data.(groups{k}); % put data into d
    Nsubj = length(d)-1; % number of subjects
    freqsX = d{1}.(block_name{1}).freqX; % x frequencies
    freqsY = d{1}.(block_name{1}).freqY; % y frequencies
    ampX = d{1}.(block_name{1}).ampX; % x amplitudes
    ampY = d{1}.(block_name{1}).ampY; % y amplitudes
    xAxis = d{1}.(block_name{1}).x_axis; % x-axis frequencies for plotting
    
    % indices of frequencies of interest
    ix = d{1}.(block_name{1}).(output).phasors.x_x.index;
    iy = d{1}.(block_name{1}).(output).phasors.y_y.index;
    
    g1 = NaN(length(xAxis),length(gblocks),Nsubj);
    g2 = NaN(length(xAxis),1000,length(gblocks));
    amps.(groups{k}).x = g1; % preallocate data structure
    amps.(groups{k}).y = g1;
    a = amps.(groups{k});
    a.bootx = g2;
    a.booty = g2;
    
    % put amplitude spectra into new data structures
    for i = 1:Nsubj
        for j = 1:length(gblocks)
            dat = d{i}.(block_name{gblocks(j)});
            pX = dat.(output).x_fft;
            pY = dat.(output).y_fft;
            
            % average across trials
            a.x_all(:,j,i) = squeeze(abs(mean(pX,2))); 
            a.y_all(:,j,i) = squeeze(abs(mean(pY,2)));
        end
    end
    
    % average across subjects
    a.x = mean(a.x_all,3); 
    a.y = mean(a.y_all,3);
    
    % permute matrices to make them easier to plot
    a.x_all = permute(a.x_all,[1 3 2]); 
    a.y_all = permute(a.y_all,[1 3 2]);
    
    % generate Figures 4A and S1A
    for i = 1:length(gblocks)
        figure(i+6) 
        subplot(2,2,k); hold on
        
        % plot circles to denote target frequencies and amplitudes
        plot(freqsX,ampX,'ko','LineWidth',lw,'MarkerFaceColor',[0 0 0])
        plot(freqsY,ampY,'ko','LineWidth',lw)
        
        % plot individual subjects
%         plot(xAxis,a.x_all(:,:,i),'Color',[0 0 0 0.3],'LineWidth',0.25...
%             ,'HandleVisibility','off') 

        % plot average avross subjects
        plot(xAxis,a.x(:,i),'Color',col(1,:),'LineWidth',lw) 
        
        % plot value of spectra at certain frequencies
        plot(freqsX,a.x(ix,i),'-ok','LineWidth',lw,'MarkerFaceColor'...
            ,[0 0 0]) 
        plot(freqsY,a.x(iy,i),'-ok','LineWidth',lw)
        
        set(gca,'Xtick',0:1:2,'Ytick',0:0.01:0.03,'box','off'...
            ,'LineWidth',1,'TickDir','out','XMinorTick','on')
        ax = gca;
        ax.XAxis.MinorTickValues = 0:0.25:2.25;
        axis([0 2.3 0 ampX(1)+.01])
        ylabel('X Amplitude (m)')
        title([graph_name{gblocks(i)}, names{k}])
        if k == 2
            legend({'X target','Y target','Hand'})
            legend boxoff
        end
        
        subplot(2,2,k+2); hold on
        
        % plot circles to denote target frequencies and amplitudes
        plot(freqsX, ampX,'ko','LineWidth',lw,'MarkerFaceColor',[0 0 0])
        plot(freqsY, ampY,'ko','LineWidth',lw)
        
        % plot individual subjects
%         plot(xAxis,a.y_all(:,:,i),'Color',[0 0 0 0.3],'LineWidth',0.25)

        % plot average avross subjects
        plot(xAxis,a.y(:,i),'LineWidth',lw,'Color',col(1,:))
        
        % plot value of spectra at certain frequencies
        plot(freqsY,a.y(iy,i),'-ok','LineWidth',lw)
        plot(freqsX,a.y(ix,i),'-ok','LineWidth',lw,'MarkerFaceColor'...
            ,[0 0 0])
        
        set(gca,'Xtick',0:1:2,'Ytick',0:0.01:0.03,'box','off'...
            ,'LineWidth',1,'TickDir','out','XMinorTick','on')
        ax = gca;
        ax.XAxis.MinorTickValues = 0:0.25:2.25;
        axis([0 2.3 0 ampY(1)+.01])
        ylabel('Y Amplitude (m)')
        xlabel('Frequency (Hz)')
    end
    
    % generate Figure S2A
    if k == 1
        subj = 4;
    else
        subj = 9;
    end
    for i = 1:length(gblocks)
        figure(i+10) % generate Figures 4A and S1A
        subplot(2,2,k); hold on
        plot(freqsX,ampX,'ko','LineWidth',lw,'MarkerFaceColor',[0 0 0])
        plot(freqsY,ampY,'ko','LineWidth',lw)
        plot(xAxis,a.x_all(:,subj,i),'Color',col(1,:),'LineWidth',lw)
        plot(freqsX,a.x_all(ix,subj,i),'-ok','LineWidth',lw...
            ,'MarkerFaceColor',[0 0 0])
        plot(freqsY,a.x_all(iy,subj,i),'-ok','LineWidth',lw)
        set(gca,'Xtick',0:1:2,'Ytick',0:0.01:0.03,'box','off'...
            ,'LineWidth',1,'TickDir','out','XMinorTick','on')
        ax = gca;
        ax.XAxis.MinorTickValues = 0:0.25:2.25;
        axis([0 2.3 0 ampX(1)+.01])
        ylabel('X Amplitude (m)')
        title([graph_name{gblocks(i)}, names{k}])
        if k == 2
            legend({'X target','Y target','Hand'})
            legend boxoff
        end
        
        subplot(2,2,k+2); hold on
        plot(freqsX, ampX,'ko','LineWidth',lw,'MarkerFaceColor',[0 0 0])
        plot(freqsY, ampY,'ko','LineWidth',lw)
        plot(xAxis,a.y_all(:,subj,i),'LineWidth',lw,'Color',col(1,:))
        plot(freqsY,a.y_all(iy,subj,i),'-ok','LineWidth',lw)
        plot(freqsX,a.y_all(ix,subj,i),'-ok','LineWidth',lw...
            ,'MarkerFaceColor',[0 0 0])
        set(gca,'Xtick',0:1:2,'Ytick',0:0.01:0.03,'box','off'...
            ,'LineWidth',1,'TickDir','out','XMinorTick','on')
        ax = gca;
        ax.XAxis.MinorTickValues = 0:0.25:2.25;
        axis([0 2.3 0 ampY(1)+.01])
        ylabel('Y Amplitude (m)')
        xlabel('Frequency (Hz)')
    end
end
end