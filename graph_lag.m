

ax = {'x_x','y_y'};
block = {'baseline','pert1','pert2','pert3','pert4','post'};
group = {'rot','mir'};
groupNames = {'Rotation','Mirror Reversal'};
Nsubj = 10;
col = lines;
col = [34 181 115
       198 156 109]./255;
% col = copper;
% col = col(floor((size(col,1)/(7))*(1:7)),:);


% store complex ratios in "ratio"
clear ratio
for i = 1:length(group)
    for j = 1:Nsubj
        for k = 1:length(block)
            for m = 1:length(ax)
                ratio(:,:,m,k,j,i) = data1.(group{i}){j}.(block{k}).cursor.phasors.(ax{m}).ratio;
            end
        end
    end
end

% average of ratio across trials
ratioMu = mean(ratio,2);

% compute phase of ratios and perform unwrapping
phase = squeeze(angle(ratioMu));
phase = unwrap(phase,[],1);
phase(5:7,2,5,3,2) = phase(5:7,2,5,3,2) - 2*pi; % correct unwrapping errors
phase(7,2,5,3,2) = phase(7,2,5,3,2) - 2*pi;

% order x- and y-phases by frequency
phase = reshape(permute(phase,[2 1 3 4 5]),[14 6 10 2]);

% calculate factor to divide phase by to get lag in seconds
freqs = sort([data1.rot{1}.baseline.freqX'; data1.rot{1}.baseline.freqY']);
scale = freqs*2*pi;
scale = repmat(scale,[1 6 10 2]);

% compute lag, convert to milliseconds, and subtract 100 ms for system
% delay
lag = -1000*phase./scale - 100;

% % order x- and y-lags by frequency
% lag = reshape(permute(lag,[2 1 3 4 5]),[14 6 10 2]);

% average across participants
lagMu = squeeze(mean(lag,3));

% plot absolute lag at baseline and late learning
col = lines;
figure(1); clf
for i = 1:2
    subplot(1,2,i); hold on
    plot(freqs,squeeze(lag(:,1,:,i)),'Color',[col(1,:) 0.5],'HandleVisibility','off')
    plot(freqs,lagMu(:,1,i),'Color',col(1,:),'LineWidth',2)
    plot(freqs,squeeze(lag(:,5,:,i)),'Color',[col(3,:) 0.5],'HandleVisibility','off')
    plot(freqs,lagMu(:,5,i),'Color',col(3,:),'LineWidth',2)
    title(groupNames{i})
    axis([0 2.2 0 1400])
    set(gca,'TickDir','out')
    xlabel('Frequency (Hz)')
    if i == 1
        ylabel('Phase lag (ms)')
    else
        legend({'Baseline','Late'})
    end
end

lagDiff = squeeze(lag(:,5,:,:) - lag(:,1,:,:));
lagDiffMu = squeeze(mean(lagDiff,2));

groupNames = {'Rotation','Mirror Reversal'};
col = [0 0 0
       237 30 121]./255;

% plot difference in lag between baseline and late learning
figure(2); clf; hold on
for i = 1:2
    plot(freqs,lagDiff(:,:,i),'Color',[col(i,:) 0.5],'HandleVisibility','off')
    plot(freqs,lagDiffMu(:,i),'Color',col(i,:),'LineWidth',2)
    axis([0 2.2 -300 1000])
    yticks(-200:200:1000)
    xticks(0:2)
    xlabel('Frequency (Hz)')
    if i == 1
        ylabel('Difference in phase lag (ms)')
    end
end
legend({'Rotation','Mirror Reversal'})
set(gca,'TickDir','out')