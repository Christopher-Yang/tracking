% plots the right hand's amplitude spectra in Figure 4A and Figure
% 4-supplements 1A and 2

function graph_ampSpectra(data)

% set variables for plotting
output = 'Rhand';
groups = {'rot','mir'};
block_name = {'baseline','pert1','pert2','pert3','pert4','post'};
graph_name = {'Baseline','Early','Train2','Train3','Late','Post'};
gblocks = [1 2 5 6]; % blocks to plot (baseline, early, late, post)
trialIdx = [1 1 8 1]; % trials within each block to plot (1 = first trial; 8 = last trial)
col = [34 181 115
       198 156 109]./255; % colors used for plotting
names = {' (Rotation)',' (Mirror Reversal)'};

for k = 1:length(groups)
    % set more variables for plotting
    d = data.(groups{k}); % put data into d
    Nsubj = length(d); % number of subjects
    freqsX = d{1}.(block_name{1}).freqX; % x frequencies
    freqsY = d{1}.(block_name{1}).freqY; % y frequencies
    ampX = d{1}.(block_name{1}).ampX; % x amplitudes
    ampY = d{1}.(block_name{1}).ampY; % y amplitudes
    xAxis = d{1}.(block_name{1}).x_axis; % x-axis frequencies for plotting
    
    % indices of frequencies of interest
    ix = d{1}.(block_name{1}).(output).phasors.x_x.index;
    iy = d{1}.(block_name{1}).(output).phasors.y_y.index;
    
    % preallocate variables
    g1 = NaN(length(xAxis),length(gblocks),Nsubj);
    g2 = NaN(length(xAxis),1000,length(gblocks));
    amps.(groups{k}).x = g1;
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
            
            % take modulus of ffts
            a.x_all(:,j,i) = squeeze(abs(pX(:,trialIdx(j)))); 
            a.y_all(:,j,i) = squeeze(abs(pY(:,trialIdx(j))));
        end
    end
    
    % average across subjects
    a.x = mean(a.x_all,3); 
    a.y = mean(a.y_all,3);
    
    % permute matrices to make them easier to plot
    a.x_all = permute(a.x_all,[1 3 2]); 
    a.y_all = permute(a.y_all,[1 3 2]);
    
    for i = 1:length(gblocks) % loop for each block to plot
        for j = 1:2 % loop for making two different figures
            if j == 1 % plots Figure 4A and Figure 4-supplement 1
                fig = [9 10];
                datX = a.x;
                datY = a.y;
            else % plots Figure 4-supplement 2
                fig = [11 12];
                if k == 1
                    subj = 4;
                else
                    subj = 2;
                end
                datX = squeeze(a.x_all(:,subj,:));
                datY = squeeze(a.y_all(:,subj,:));
            end
                
            % plot x-axis spectra
            figure(fig(1))
            subplot(4,2,2*(i-1)+k); hold on
            
            % plot target frequencies and amplitudes
            plot(freqsX,ampX*100,'d','Color',col(1,:),'MarkerSize',6,...
                'MarkerFaceColor',col(1,:))
            plot(freqsY,ampY*100,'d','Color',col(2,:),'MarkerSize',6,...
                'MarkerFaceColor',col(2,:))
            
            % plot value of spectra at certain frequencies
            plot(freqsX,datX(ix,i)*100,'-o','Color',col(1,:),...
                'MarkerSize',6,'MarkerFaceColor',col(1,:))
            plot(freqsY,datX(iy,i)*100,'-o','Color',col(2,:),...
                'MarkerSize',6,'MarkerFaceColor',col(2,:))
            
            % plot average across subjects
            plot(xAxis,datX(:,i)*100,'k','LineWidth',0.5)
            
            set(gca,'Xtick',0:2,'Ytick',0:3,'box','off','TickDir','out')
            set(gcf,'Position',[100 50 520 600])
            axis([0 2.3 0 2.5])
            title([graph_name{gblocks(i)}, names{k}])
            ylabel('X Amplitude (cm)')
            if i == length(gblocks)
                xlabel('Frequency (Hz)')
            elseif i == 1 && k == 2
                legend({'X target','Y target','Hand response (X freq)',...
                    'Hand response (Y freq)'})
                legend boxoff
            end
            
            % plot y-axis spectra
            figure(fig(2))
            subplot(4,2,2*(i-1)+k); hold on
            
            % plot target frequencies and amplitudes
            plot(freqsX, ampX*100,'d','Color',col(1,:),'MarkerSize',6,...
                'MarkerFaceColor',col(1,:))
            plot(freqsY, ampY*100,'d','Color',col(2,:),'MarkerSize',6,...
                'MarkerFaceColor',col(2,:))
            
            % plot value of spectra at certain frequencies
            plot(freqsY,datY(iy,i)*100,'-o','Color',col(2,:),...
                'MarkerSize',6,'MarkerFaceColor',col(2,:))
            plot(freqsX,datY(ix,i)*100,'-o','Color',col(1,:),...
                'MarkerSize',6,'MarkerFaceColor',col(1,:))
            
            % plot individual subjects
%             plot(xAxis,a.y_all(:,:,i),'Color',[0 0 0 0.3],...
%                 'LineWidth',0.25)
            
            % plot average avross subjects
            plot(xAxis,datY(:,i)*100,'k','LineWidth',0.5)
            
            set(gca,'Xtick',0:2,'Ytick',0:3,'box','off','TickDir','out')
            set(gcf,'Position',[100 50 520 600])
            axis([0 2.3 0 2.5])
            title([graph_name{gblocks(i)}, names{k}])
            ylabel('Y Amplitude (cm)')
            if i == length(gblocks)
                xlabel('Frequency (Hz)')
            elseif i == 1 && k == 2
                legend({'X target','Y target','Hand response (X freq)',...
                    'Hand response (Y freq)'})
                legend boxoff
            end
        end
    end
end
end