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
Ntrial = length(data.day2{1}.B1_baseline.MSE);

% colors for plotting
col = [180 180 0
       0 191 255
       255 99 71
       251 176 59
       140 98 57
       198 156 109]./255;

% store all phasors
for i = 1:Ngroup
    for j = 1:allSubj(i)
        for k = 1:Nblock
            phasor{i}(:,j,k) = mean(data.(group{i}){j}.(block{i}{k}).cursor.phasors.x_x.ratio,2); % average across trials
            phasor_all{i}(:,:,j,k) = data.(group{i}){j}.(block{i}{k}).cursor.phasors.x_x.ratio; % no averaging across trials
        end
    end
end

% compute normalized gain for individual subjects
percent = NaN(Nfreq,Ntrial,max(allSubj),Ngroup,2);
for m = 1:2 % loop over two flip blocks
    
    for k = 1:Ngroup % loop over groups
        Nsubj = allSubj(k); % number of subjects in current group
        
        for j = 1:Nsubj % loop over subjects
            
            for i = 1:Nfreq % loop over frequencies
                
                for n = 1:Ntrial
                    mu = [real(phasor_all{k}(i,n,j,1)) imag(phasor_all{k}(i,n,j,1))]; % late learning phasor
                    mu2 = [real(phasor_all{k}(i,n,j,m+1)) imag(phasor_all{k}(i,n,j,m+1))]; % flip block phasor
                    
                    unitVec = mu/norm(mu); % unit vector of late learning phasor
                    flip = dot(mu2,unitVec); % gain of flip block phasor projected onto late learning phasor
                    late = norm(mu); % gain of late learning pahsor
                    percent(i,n,j,k,m) = flip/late; % normalize gain based on late learning
                end
            end
        end
    end
end

% average data from highest three frequencies
highFreq = percent(4:6,:,:,:,:);
flip = {'flip1','flip2'};
for k = 1:Ngroup
    for j = 1:2
        gain.(group{k}).(flip{j}) = squeeze(mean(mean(highFreq(:,:,1:allSubj(k),k,j),1),2));
    end
end

% load habit weights from mixture model for point-to-point task
load Variables/weight2_opt

%% Supplementary Fig 3A

% index of subjects to be plotted
idx{1} = [2 6 7 10 11];
idx{2} = [5 7 12 13 14];
idx{3} = 1:5;

f = figure(7); clf
set(f,'Position',[200 200 500 600]);
for k = 1:Ngroup
    for j = 1:5
        subplot(3,5,(k-1)*5+j); hold on
        plot([0 2],[0 0],'k')
        plot([0 2],[1 1],'--','Color',col(4,:))
        plot([0 2],[-1 -1],'--','Color',col(5,:))
        plot(freq,percent(:,:,idx{k}(j),k,1),'.','Color',col(k,:),'MarkerSize',8)
        plot(freq,mean(percent(:,:,idx{k}(j),k,1),2),'-ko','MarkerFaceColor',col(k,:),'MarkerSize',4,'LineWidth',1)
        axis([0 1.65 -3 2])
        set(gca,'TickDir','out')
        if k == 2 && j == 1
            ylabel('Gain')
        elseif k == 3 && j == 3
            xlabel('Frequency (Hz)')
        end
    end
end

% save figure for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/gain_single_subj','-dpdf','-painters')

%% Figure 5C
% linear regression between tracking and p2p tasks, excluding 1 subject 
% from the 10-day group
x = [gain.day2.flip1; gain.day5.flip1; gain.day10.flip1(2:5)];
y = [weight2_opt{1}(:,4); weight2_opt{2}(:,4); weight2_opt{3}(2:5,4)];
p = polyfit(x,y,1);
yfit = polyval(p,[-1.3; x]);
lm = fitlm(x,y);

% linear regression between tracking and p2p tasks with all subjects
x2 = [gain.day2.flip1; gain.day5.flip1; gain.day10.flip1];
y2 = [weight2_opt{1}(:,4); weight2_opt{2}(:,4); weight2_opt{3}(:,4)];
p2 = polyfit(x2,y2,1);
yfit2 = polyval(p2,[-1.3; x2]);
lm2 = fitlm(x2,y2);

% display results of regression
disp('Stats for Fig 5C: linear regression between tracking and p2p tasks')
disp(['   slope = ' num2str(lm.Coefficients.Estimate(2))])
disp(['   p = ' num2str(lm.Coefficients.pValue(2))])
disp(['   r = ' num2str(sqrt(lm.Rsquared.Ordinary)) ])

disp('Stats for linear regression between tracking and p2p tasks no outliers excluded')
disp(['   slope = ' num2str(lm2.Coefficients.Estimate(2))])
disp(['   p = ' num2str(lm2.Coefficients.pValue(2))])
disp(['   r = ' num2str(sqrt(lm2.Rsquared.Ordinary))])

% plot figure
f = figure(8); clf; hold on
set(f,'Position',[200 200 160 120]);
plot([-1.3; x],yfit,'k','LineWidth',1,'HandleVisibility','off')
plot(gain.day2.flip1,weight2_opt{1}(:,4),'.','Color',col(1,:),'MarkerSize',10)
plot(gain.day5.flip1,weight2_opt{2}(:,4),'.','Color',col(2,:),'MarkerSize',10)
plot(gain.day10.flip1,weight2_opt{3}(:,4),'.','Color',col(3,:),'MarkerSize',10)
plot(gain.day10.flip1(1),weight2_opt{3}(1,4),'xk','LineWidth',0.5,'MarkerSize',10)
xlabel('Habitual gain (tracking)')
ylabel('Probability of habitual reach (point-to-point)')
xticks(-2:0)
yticks(0:0.25:0.75)
axis([-2 0.5 0 0.75])

% save figure for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/scatter','-dpdf','-painters')

%% statistical analysis for single subjects
% variables ending in 1 correspond to first flip block
% variables ending in 2 correspond to second flip block
pvalues1 = NaN(6,1);
pvalues2 = NaN(6,1);
habit1 = NaN(max(allSubj),3);
habit2 = NaN(max(allSubj),3);
for k = 1:Ngroup
    for j = 1:allSubj(k)
        flag1 = 0; % flags for identifying significant comparisons
        flag2 = 0;
        for i = 1:Nfreq
            % one-sample t-tests
            [~, p] = ttest(percent(i,:,j,k,1),0,'Tail','left');
            pvalues1(i) = p;
            [~, p] = ttest(percent(i,:,j,k,2),0,'Tail','left');
            pvalues2(i) = p;
        end
        
        % apply Holm-Bonferroni correction
        pvalues1 = sort(pvalues1);
        pvalues2 = sort(pvalues2);
        for i = 1:Nfreq
            if flag == 1
                continue
            end
            
            % check to see if p-value is under a corrected alpha
            if pvalues1(i) < 0.05/(Nfreq + 1 - i)
                flag1 = 1;
            end
            if pvalues2(i) < 0.05/(Nfreq + 1 - i)
                flag2 = 1;
            end
        end
        
        % store results of statistical comparisons
        habit1(j,k) = flag1;
        habit2(j,k) = flag2;
    end
end

% remove data from subject 1 in 10-day group
habit1(1,3) = NaN;
habit2(1,3) = NaN;

% number of subjects in each group that exhibited habitual behavior
Nhabit1 = sum(habit1,'omitnan');
Nhabit2 = sum(habit2,'omitnan');

disp('Stats for Fig 5B: # subjects with significantly negative gains')
disp('(one-sample t-test with Holm-Bonferroni correction)')
disp('   Flip 1')
disp(['      2-day: ' num2str(Nhabit1(1)) ' of 13'])
disp(['      5-day: ' num2str(Nhabit1(2)) ' of 14'])
% disp(['      10-day: ' num2str(Nhabit1(3)) ' of 4'])
disp('   Flip 2')
disp(['      2-day: ' num2str(Nhabit2(1)) ' of 13'])
disp(['      5-day: ' num2str(Nhabit2(2)) ' of 14'])
% disp(['      10-day: ' num2str(Nhabit2(3)) ' of 4'])

%% group level analysis
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

% store data from outlier subject
outlier = squeeze(percent(:,1,3,:));

% remove outlier (subject 1 from 10-day group)
percent(:,1,3,:) = NaN;

% mean and standard error of gains
percent_mu = squeeze(mean(percent,2,'omitnan'));
percent_se = squeeze(std(percent,[],2,'omitnan'))./sqrt(repmat([13 14 4],[6 1 2]));

% save data for analysis in R
y = percent(~isnan(percent(:,:,1:2,1)));
groupNames(1:78,1) = "2-day";
groupNames(79:162,1) = "5-day";
frequency = repmat((1:Nfreq)',[sum([13 14]) 1]);
s1 = repmat(1:13,[Nfreq 1]);
s2 = repmat(14:27,[Nfreq 1]);
subject = [s1(:); s2(:)];
T = table(groupNames, frequency, subject, y, 'VariableNames', {'group','frequency','subject','gain'});
% writetable(T,'C:/Users/Chris/Documents/R/habit/data/habitGain.csv')

% statistical analysis for group-level data
pvalues1 = NaN(6,1);
pvalues2 = NaN(6,1);
habit1 = zeros(Nfreq,Ngroup);
habit2 = zeros(Nfreq,Ngroup);
for k = 1:Ngroup
    
    % one-sample t-tests
    for i = 1:Nfreq
        [~, p] = ttest(percent(i,1:allSubj(k),k,1),0,'Tail','left');
        pvalues1(i) = p;
        [~, p] = ttest(percent(i,1:allSubj(k),k,2),0,'Tail','left');
        pvalues2(i) = p;
    end
    
    % apply Holm-Bonferroni correction
    [pvalues1, idx] = sort(pvalues1);
    [pvalues2, idx2] = sort(pvalues2);
    for i = 1:Nfreq
        if pvalues1(i) < 0.05/(Nfreq + 1 - i)
            habit1(i,k) = 1;
        end
        if pvalues2(i) < 0.05/(Nfreq + 1 - i)
            habit2(i,k) = 1;
        end
    end
    
    % order of frequencies analyzed for habit1 and habit2
    freqIdx(:,k) = idx;
    freqIdx2(:,k) = idx2;
end

Nhabit1 = sum(habit1);
Nhabit2 = sum(habit2);

disp('Stats for Fig 5B: # frequencies with significantly negative gains')
disp('(one-sample t-test with Holm-Bonferroni correction)')
disp('   Flip 1')
disp(['      2-day: ' num2str(Nhabit1(1)) ' of 6'])
disp(['      5-day: ' num2str(Nhabit1(2)) ' of 6'])
% disp(['      10-day: ' num2str(Nhabit1(3)) ' of 6'])
disp('   Flip 2')
disp(['      2-day: ' num2str(Nhabit2(1)) ' of 6'])
disp(['      5-day: ' num2str(Nhabit2(2)) ' of 6'])
% disp(['      10-day: ' num2str(Nhabit2(3)) ' of 6'])

%% Supplementary Figure 3B

f = figure(10); clf
set(f,'Position',[200 200 350 170]);
for i = 1:2
    subplot(1,2,i); hold on
    plot([0 2],[0 0],'k','HandleVisibility','off')
    plot([0 2],[1 1],'--','Color',col(4,:),'HandleVisibility','off')
    plot([0 2],[-1 -1],'--','Color',col(5,:),'HandleVisibility','off')
    
    for j = 1:Ngroup
        plot(freq, percent(:,:,j,i), 'Color', [col(j,:) 0.5])
        plot(freq, percent_mu(:,j,i), '-o', 'Color', col(j,:), 'MarkerFaceColor', col(j,:), 'LineWidth', 1.5, 'MarkerSize', 5)
    end
    
    plot(freq, outlier(:,i), '--', 'Color', col(j,:))
    
    set(gca,'TickDir','out')
    axis([0 1.6 -2.1 1.2])
    xlabel('Frequency (Hz)')
    xticks(0:0.5:2)
    yticks(-2:1)
    if i == 1
        title('Flip 1')
        ylabel('Gain')
    else
        title('Flip 2')
    end
end
% save figure for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/gain_habit3','-dpdf','-painters')

%% Figure 5B
f = figure(9); clf
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

end