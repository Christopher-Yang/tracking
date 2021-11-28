function graph_flipGain(data)

addpath ../shadedErrorBar

group = {'day2','day5','day10'};
names = {'2-day','5-day','10-day'};
allSubj = [13 14 5];
Ngroup = 3;
block{1} = {'B9', 'B11_habit','B12_habit'};
block{2} = {'B24', 'B26_habit','B27_habit'};
block{3} = {'B49', 'B51_habit','B52_habit'};
% block{1} = {'B10_dark', 'B13_habitDark'};
% block{2} = {'B25_dark', 'B28_habitDark'};
% block{3} = {'B50_dark', 'B53_habitDark'};
freq = data.day2{1}.B1_baseline.freqX;

Nfreq = length(freq);
col = copper;
col = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

Nblock = length(block{1});

for i = 1:Ngroup
    for j = 1:allSubj(i)
        for k = 1:Nblock
            phasor{i}(:,j,k) = mean(data.(group{i}){j}.(block{i}{k}).cursor.phasors.x_x.ratio,2);
        end
    end
end

%%
figure(1); clf
for k = 1:Nblock
    for j = 1:Ngroup
        subplot(Nblock,Ngroup,3*(k-1) + j); hold on
        plot([0 0], [-1 1], 'k')
        plot([-1 1], [0 0], 'k')
        for i = 1:6
            plot(real(phasor{j}(i,:,k)), imag(phasor{j}(i,:,k)), '.', 'MarkerSize', 10, 'Color', col(i,:))
            plot(real(mean(phasor{j}(i,:,k), 2)), imag(mean(phasor{j}(i,:,k), 2)), '.', 'MarkerSize', 30, 'Color', col(i,:))
        end
        axis([-1 1 -1 1])
        axis square
        if k == 1
            title(names{j})
            if j == 1
                ylabel('Late')
            end
        elseif k == 2 && j == 1
            ylabel('Flip 1')
        elseif k == 3 && j == 1
            ylabel('Flip 2')
        end
    end
end

%%
figure(2); clf
for j = 1:3
    for i = 1:6
        mu = [mean(real(phasor{j}(i,:,1))) mean(imag(phasor{j}(i,:,1)))];
        mu2 = [mean(real(phasor{j}(i,:,2))) mean(imag(phasor{j}(i,:,2)))];
        
        slope = mu(2)/mu(1);
        y = [-1*slope slope];
        unitVec = mu/norm(mu);
        len = dot(mu2,unitVec);
        
        project = len*[unitVec(1) unitVec(2)];
        
        subplot(3,6,(j-1)*6 + i); hold on
        plot([0 0], [-1 1], 'k')
        plot([-1 1], [0 0], 'k')
        plot(real(phasor{j}(i,:,1)),imag(phasor{j}(i,:,1)),'.k','MarkerSize',10)
        plot(real(phasor{j}(i,:,2)),imag(phasor{j}(i,:,2)),'.r','MarkerSize',10)
        plot(mu(1),mu(2),'.k','MarkerSize',30)
        plot(mu2(1),mu2(2),'.r','MarkerSize',30)
        plot([mu2(1) project(1)],[mu2(2) project(2)],'r')
        plot([-1 1],y,'k')
        axis([-1 1 -1 1])
        axis square
    end
end

%%

percent = NaN(Nfreq,max(allSubj),Ngroup,2);
for m = 1:2
    for k = 1:Ngroup
        Nsubj = allSubj(k);
        for j = 1:Nsubj
            for i = 1:Nfreq
                
                mu = [real(phasor{k}(i,j,1)) imag(phasor{k}(i,j,1))];
                mu2 = [real(phasor{k}(i,j,m+1)) imag(phasor{k}(i,j,m+1))];
                
                unitVec = mu/norm(mu);
                flip = dot(mu2,unitVec);
                
                late = norm(mu);
                
                percent(i,j,k,m) = flip/late;
            end
        end
    end
end

percent_mu = squeeze(mean(percent,2,'omitnan'));
percent_se = squeeze(std(percent,[],2,'omitnan'))./sqrt(repmat(allSubj,[6 1 2]));

y = percent(~isnan(percent));

groupNames([1:78 193:270],1) = "2-day";
groupNames([79:162 271:354],1) = "5-day";
groupNames([163:192 355:384],1) = "10-day";
frequency = repmat((1:Nfreq)',[sum(allSubj)*2 1]);
blockNames(1:192,1) = "Flip1";
blockNames(193:384,1) = "Flip2";
s1 = repmat(1:13,[Nfreq 1]);
s2 = repmat(14:27,[Nfreq 1]);
s3 = repmat(28:32,[Nfreq 1]);
subject = [s1(:); s2(:); s3(:); s1(:); s2(:); s3(:)];
T = table(groupNames, frequency, blockNames, subject, y, 'VariableNames', {'group','frequency','block', 'subject','gain'});
writetable(T,'C:/Users/Chris/Documents/R/habit/data/habitGain.csv')

%%

% col1 = [0 60 0]/255;
% col2 = [50 255 50]/255;
% map1 = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2),Nfreq)', linspace(col1(3),col2(3),Nfreq)'];

% col1 = [128 0 128]/255;
% col2 = [230 230 250]/255;
% map2 = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2),Nfreq)', linspace(col1(3),col2(3),Nfreq)'];

labels = {'Flip 1','Flip 2'};
f = figure(3); clf
set(f,'Position',[200 200 250 120]);
for j = 1:Ngroup
    subplot(1,3,j); hold on
    plot([0 3],[1 1],'r--','HandleVisibility','off')
    plot([0 3],[0 0],'k','HandleVisibility','off')
    plot([0 3],[-1 -1],'g--','HandleVisibility','off')
    for i = 1:Nfreq
        s = shadedErrorBar(1:2,squeeze(percent_mu(i,j,:)),squeeze(percent_se(i,j,:)),'lineprops','-o');
        editErrorBar(s,col(i,:),1);
    end
    axis([0.8 2.2 -1.2 1.2]) 
    xticks(1:2)
    xticklabels(labels)
    yticks(-1:1)
    if j == 1
        ylabel('Gain relative to late learning')
    end
    set(gca,'TickDir','out')
    title(names{j})
end

% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/gain_habit','-dpdf','-painters')

%%

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
        ylabel('Gain relative to late learning')
    else
        title('Flip 2')
    end
end
print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/gain_habit2','-dpdf','-painters')


end