function graph_lag(data1)
% plots the increase in lag between hand and target movements from baseline
% to late learning

% set variables for analysis
ax = {'x_x','y_y'}; % axes in which to evaluate lag; first letter is input axis, second letter is output axis
block = {'baseline','pert1','pert2','pert3','pert4','post'}; % blocks to analyze 
group = {'rot','mir'}; % names of groups
groupNames = {'Rotation','Mirror Reversal'}; % group names for creating figure legends
Nsubj = 10; % number of subjects

% colors for plotting data
col = lines;
col2 = [0 0 0
       237 30 121]./255;

% store complex ratios in "ratio"
clear ratio
for i = 1:length(group)
    for j = 1:Nsubj
        for k = 1:length(block)
            for m = 1:length(ax)
                ratio(:,:,m,k,j,i) = data1.(group{i}){j}.(block{k}). ...
                    cursor.phasors.(ax{m}).ratio;
            end
        end
    end
end

% average of ratio across trials
ratioMu = mean(ratio,2);

% compute phase of ratios and perform unwrapping
phase = squeeze(angle(ratioMu));
phase = unwrap(phase,[],1);
phase(5:7,2,5,3,2) = phase(5:7,2,5,3,2) - 2*pi; % corrects unwrapping errors
phase(7,2,5,3,2) = phase(7,2,5,3,2) - 2*pi; % corrects unwrapping errors

% order x- and y-phases by frequency
phase = reshape(permute(phase,[2 1 3 4 5]),[14 6 10 2]);

% calculate factor to divide phase by to get lag in seconds
freqs = sort([data1.rot{1}.baseline.freqX'; data1.rot{1}.baseline.freqY']);
scale = freqs*2*pi;
scale = repmat(scale,[1 6 10 2]);

% compute lag, convert to milliseconds, and subtract 100 ms for system
% delay
lag = -1000*phase./scale - 100;

% average across participants
lagMu = squeeze(mean(lag,3));

% compute increase in lag from baseline to late learning
lagDiff = squeeze(lag(:,5,:,:) - lag(:,1,:,:));
lagDiffMu = squeeze(mean(lagDiff,2));

% plot figures
figure(14); clf
for i = 1:2
    % absolute lag at baseline and late learning; Figure 4-figure 
    % supplement 1C
    subplot(2,2,i); hold on
    plot(freqs,squeeze(lag(:,1,:,i)),'Color',[col(1,:) 0.5],...
        'HandleVisibility','off')
    plot(freqs,lagMu(:,1,i),'Color',col(1,:),'LineWidth',2)
    plot(freqs,squeeze(lag(:,5,:,i)),'Color',[col(3,:) 0.5],...
        'HandleVisibility','off')
    plot(freqs,lagMu(:,5,i),'Color',col(3,:),'LineWidth',2)
    title(groupNames{i})
    axis([0 2.2 0 1400])
    set(gca,'TickDir','out')
    xlabel('Frequency (Hz)')
    if i == 1
        ylabel('Lag between hand and target(ms)')
    else
        legend({'Baseline','Late'})
    end

    % increase in lag between baseline and late learning; Figure 4C
    subplot(2,2,3:4); hold on
    plot(freqs,lagDiff(:,:,i),'Color',[col2(i,:) 0.5],...
        'HandleVisibility','off')
    plot(freqs,lagDiffMu(:,i),'Color',col2(i,:),'LineWidth',2)
    if i == 1
        axis([0 2.2 -300 1000])
        yticks(-200:200:1000)
        xticks(0:2)
        xlabel('Frequency (Hz)')
        title('Late learning lag - baseline lag')
        ylabel('Increase in lag (ms)')
    end
end
legend({'Rotation','Mirror Reversal'})
set(gca,'TickDir','out')

end