% graph_flipGain plots the gain of cursor movements from the flip block
% normalized to the gain at late learning

function graph_flipGain(data)

% set variables for analysis
group = {'day2','day5','day10'}; % name of groups
allSubj = [length(data.day2) length(data.day5) length(data.day10)]; % number of subjects in each group
block{1} = {'B9', 'B11_habit','B12_habit'}; % block names
block{2} = {'B24', 'B26_habit','B27_habit'};
block{3} = {'B49', 'B51_habit','B52_habit'};
freq = data.day2{1}.B1_baseline.freqX; % x-axis target frequencies
Ngroup = length(group); % number of groups
Nfreq = length(freq); % number of frequencies
Nblock = length(block{1}); % number of blocks

% store all phasors into "phasor"
for i = 1:Ngroup
    for j = 1:allSubj(i)
        for k = 1:Nblock
            phasor{i}(:,j,k) = mean(data.(group{i}){j}.(block{i}{k}).cursor.phasors.x_x.ratio,2);
        end
    end
end

% compute gain of flip block normalized to late learning gain
percent = NaN(Nfreq,max(allSubj),Ngroup,2); % preallocate variable
for m = 1:2 % loop over two flip blocks
    
    for k = 1:Ngroup % loop over groups
        Nsubj = allSubj(k); % number of subjects in current group
        
        for j = 1:Nsubj % loop over subjects
            
            for i = 1:Nfreq % loop over frequencies
                
                mu = [real(phasor{k}(i,j,1)) imag(phasor{k}(i,j,1))]; % late learning phasor
                mu2 = [real(phasor{k}(i,j,m+1)) imag(phasor{k}(i,j,m+1))]; % flip block phasor
                
                unitVec = mu/norm(mu); % unit vector of late learning phasor
                flip = dot(mu2,unitVec); % gain of flip block phasor projected onto late learning phasor
                late = norm(mu); % gain of late learning pahsor
                percent(i,j,k,m) = flip/late; % normalize gain based on late learning
            end
        end
    end
end

% percent(:,1,3,:) = NaN;

% mean and standard error of gains
percent_mu = squeeze(mean(percent,2,'omitnan'));
percent_se = squeeze(std(percent,[],2,'omitnan'))./sqrt(repmat([13 14 4],[6 1 2]));

% save data for analysis in R
% y = percent(~isnan(percent));
% groupNames([1:78 187:264],1) = "2-day";
% groupNames([79:162 265:348],1) = "5-day";
% groupNames([163:186 349:372],1) = "10-day";
% frequency = repmat((1:Nfreq)',[sum([13 14 4])*2 1]);
% blockNames(1:186,1) = "Flip1";
% blockNames(187:372,1) = "Flip2";
% s1 = repmat(1:13,[Nfreq 1]);
% s2 = repmat(14:27,[Nfreq 1]);
% s3 = repmat(28:31,[Nfreq 1]);
% subject = [s1(:); s2(:); s3(:); s1(:); s2(:); s3(:)];
% T = table(groupNames, frequency, blockNames, subject, y, 'VariableNames', {'group','frequency','block', 'subject','gain'});
% writetable(T,'C:/Users/Chris/Documents/R/habit/data/habitGain.csv')

%% plot Figure 5B

% colors for plotting
col = [180 180 0
       0 191 255
       255 99 71
       251 176 59
       140 98 57]./255;

f = figure(4); clf
set(f,'Position',[200 200 250 140]);
for i = 1:2
    subplot(1,2,i); hold on
    plot([0 2],[0 0],'k','HandleVisibility','off')
    plot([0 2],[1 1],'--','Color',col(4,:),'HandleVisibility','off')
    plot([0 2],[-1 -1],'--','Color',col(5,:),'HandleVisibility','off')
    for j = 1:Ngroup
        s = shadedErrorBar(freq,percent_mu(:,j,i),percent_se(:,j,i),'lineProps','-o');
        editErrorBar(s,col(j,:),1);
    end
    set(gca,'TickDir','out')
    axis([0 1.6 -1.2 1.2])
    xlabel('Frequency (Hz)')
    xticks(0:0.5:2)
    yticks(-1:1)
    if i == 1
        title('Flip 1')
        ylabel('Gain')
    else
        title('Flip 2')
    end
end

% save figure for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/gain_habit2','-dpdf','-painters')

%%
figure(5); clf
for i = 1:2
    subplot(1,2,i); hold on
    plot([0 2],[0 0],'k','HandleVisibility','off')
    plot([0 2],[1 1],'--','Color',col(4,:),'HandleVisibility','off')
    plot([0 2],[-1 -1],'--','Color',col(5,:),'HandleVisibility','off')
    for j = 1:Ngroup
        if j == 3
            plot(freq,percent(:,2:5,j,i),'Color',[col(j,:) 0.5])
            plot(freq,mean(percent(:,2:5,j,i),2,'omitnan'),'Color',col(j,:),'LineWidth',2)
        else
            plot(freq,percent(:,:,j,i),'Color',[col(j,:) 0.5])
            plot(freq,mean(percent(:,:,j,i),2,'omitnan'),'Color',col(j,:),'LineWidth',2)
        end
    end
    axis([0 1.6 -1.2 1.2])
end


end