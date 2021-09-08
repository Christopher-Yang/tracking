function graph_gainMatrix(data, block_name, graph_name)

output = 'cursor';
groups = {'day2','day5','day10'};
groupNames = {'2-day','5-day','10-day'};
names = {'x_x','x_y','y_x','y_y'};

Ngroup = length(groups);
maxNormal = 0;
counterbalance{1} = 8:13;
counterbalance{2} = 8:14;
counterbalance{3} = 4:5;

if strcmp(output,'cursor')
    labels = {'X_{T}X_{O} (response)','X_{T}Y_{O}','Y_{T}X_{O}','Y_{T}Y_{O} (response)'};
    workspace = 1:4;
elseif strcmp(output,'Rhand')
    labels = {'X_{T}X_{O}','X_{T}Y_{O}','Y_{T}X_{O} (response)','Y_{T}Y_{O}'};
    workspace = [1 3];
    nullspace = [2 4];
    response = 3;
elseif strcmp(output,'Lhand')
    labels = {'X_{T}X_{O}','X_{T}Y_{O} (response)','Y_{T}X_{O}','Y_{T}Y_{O}'};
    workspace = [2 4];
    nullspace = [1 3];
    response = 2;
end

for q = 1:Ngroup % loop over groups
    clear opt1 opt2 template
    Nsubj = length(data.(groups{q}));
    blk = block_name.(groups{q});
    Nblock = length(blk);
    Nfreq = length(data.(groups{q}){1}.(blk{1}).freqX);
    Ntrial = 5;
    paramsInit = zeros([2 Nblock]);
    
    dark{q} = find(contains(graph_name.(groups{q}),'D)'));
    flip{q} = find(contains(graph_name.(groups{q}),'(F)'));
    special{q} = find(contains(graph_name.(groups{q}),'('));
    normal{q} = 1:Nblock;
    normal{q}(special{q}) = [];
    
    if maxNormal < length(normal{q})
        maxNormal = length(normal{q});
    end
    
%     set remove = 1 to remove special blocks
%     remove = 0;
%     if remove
%         index = normal;
%     end
    
    for p = 1:Nsubj % loop over subjects
        for k = 1:Nblock % loop over block
            for i = 1:4 % loop over all combinations of hand/target axes
                a = data.(groups{q}){p}.(blk{k}).(output).phasors.(names{i}).ratio;
                if size(a,2) ~= Ntrial
                    Nmissing = Ntrial - size(a,2);
                    cplx{q}(:,:,i,k,p) = [a NaN(Nfreq,Nmissing)];
                else
                    cplx{q}(:,:,i,k,p) = a; % put all complex ratios in cplx
                end
            end
        end
        
        for i = 1:Nfreq
            
            % compute average baseline phasor for a particular frequency
            dat = mean(cplx{q}(i,:,[1 4],1,p),2,'omitnan');

            % x-component of xx and yy baseline phasor
            x = cos(angle(dat));
            
            % y-component of xx and yy baseline phasor
            y = sin(angle(dat));
            
            % project data onto baseline x phasor
            dat = cplx{q}(i,:,[1 2],:,p); % extract data for x-target to x-/y-hand 
            phasor = reshape(dat, [1 numel(dat)]);
            num = [real(phasor); imag(phasor)]; % separate real and imaginary parts of phasors
            unitVec = [x(1); y(1)]; % unit vector with phase equal to baseline vector
            gainX{q}(:,i,p) = dot(num, repmat(unitVec, [1 numel(dat)])); % project num onto unitVec
            frac1{q}(:,i,p) = abs(gainX{q}(:,i,p))./abs(phasor)'; % calculate proportion of gain retained in projection
            
            % project data onto baseline y phasor
            dat = cplx{q}(i,:,[3 4],:,p); % extract data for x-target to x-/y-hand 
            phasor = reshape(dat, [1 numel(dat)]);
            num = [real(phasor); imag(phasor)]; % separate real and imaginary parts of phasors
            unitVec = [x(2); y(2)]; % unit vector with phase equal to baseline vector
            gainY{q}(:,i,p) = dot(num, repmat(unitVec, [1 numel(dat)])); % project num onto unitVec
            frac2{q}(:,i,p) = abs(gainY{q}(:,i,p))./abs(phasor)'; % calculate proportion of gain retained in projection
            
        end
    end
    
    % combine opt1 and opt2
    thetaOpt{q} = [reshape(gainX{q},[Ntrial 2 Nblock Nfreq Nsubj]) ...
        reshape(gainY{q},[Ntrial 2 Nblock Nfreq Nsubj])];
    thetaOpt{q} = permute(thetaOpt{q}, [2 3 4 1 5]);
    
    % need to swap signs to properly counterbalance if not looking at
    % cursor
%     thetaOpt{q}(workspace,2:end,:,:,counterbalance{q}) = -thetaOpt{q}(workspace,2:end,:,:,counterbalance{q});
    
    thetaOpt_mu{q} = mean(thetaOpt{q},5,'omitnan');
    thetaOpt_se{q} = std(thetaOpt{q},[],5,'omitnan')./sqrt(Nsubj);
    
    % combine frac1 and frac2
    fracOpt{q} = [reshape(frac1{q},[Ntrial 2 Nblock Nfreq Nsubj]) ...
        reshape(frac2{q},[Ntrial 2 Nblock Nfreq Nsubj])];
    fracOpt{q} = permute(fracOpt{q}, [2 3 4 1 5]);
    
    % weight thetaOpt by fraction of gain lost during projection
%     thetaOpt{q} = thetaOpt{q} .* (1./fracOpt{q});
    
    % shape thetaOpt into gain matrix format: the first and second dimensions
    % of gainMat correspond to the gain matrix for a given frequency (third
    % dimension), block (fourth dimension), subject (fifth dimension), and
    % group (sixth dimension)
    gainMat{q} = reshape(thetaOpt{q},[2 2 Nblock Nfreq Ntrial Nsubj]);
    gainMat{q} = permute(gainMat{q},[1 2 4 3 5 6]);
    
    habitGains = thetaOpt{q}(1,end-2,:,:,:);
    habitGains = squeeze(habitGains);
    baselineGains = thetaOpt{q}(1,1,:,:,:);
    baselineGains = squeeze(baselineGains);
    percent = habitGains./baselineGains;
    
    flipSign.(groups{q}) = percent(:,:,:,1) < -0.2;
end

save flipSign flipSign

y = [];
for k = 1:Ngroup
    for j = 1
        for i = 1:2
            if i == 1
                t = mean(thetaOpt{k}(j,1,:,:,:),4,'omitnan');
            else
                t = mean(thetaOpt{k}(j,5,:,:,:),4,'omitnan');
            end
            y = [y; t(:)];
        end
    end
end

g(1:156,1) = "2-day";
g(157:324,1) = "5-day";
g(325:384,1) = "10-day";
b([1:78 157:240 325:354],1) = "Baseline";
b([79:156 241:324 355:384],1) = "Late";
frequency = repmat((1:Nfreq)',[64 1]);
s1 = repmat(1:13,[Nfreq 2]);
s2 = repmat(14:27,[Nfreq 2]);
s3 = repmat(28:32,[Nfreq 2]);
subject = [s1(:); s2(:); s3(:)];
T = table(g, b, frequency, subject, y, 'VariableNames', {'group','block','frequency','subject','gain'});
writetable(T,'C:/Users/Chris/Documents/R/habit/data/gain_skill.csv')

% for plotting lines
col = lines;
col = col(1:7,:);

% for plotting heatmaps
col1 = [1 0 0];
col2 = [1 1 1];
Nstep = 100;
map1 = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2),Nstep)', linspace(col1(3),col2(3),Nstep)'];

col1 = [1 1 1];
col2 = [0 0 1];
map2 = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2),Nstep)', linspace(col1(3),col2(3),Nstep)'];

map = [map1; map2];
clims = [-1 1];

% for plotting vectors
col1 = [0 128 0]/255;
col2 = [152 251 152]/255;
map1 = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2),Nfreq)', linspace(col1(3),col2(3),Nfreq)'];

col1 = [128 0 128]/255;
col2 = [230 230 250]/255;
map2 = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2),Nfreq)', linspace(col1(3),col2(3),Nfreq)'];

%% plot vectors

labels = {'Baseline','Early','Late',"Flip 1","Flip 2"};
figure(1); clf
for q = 1:Ngroup
    blk = [1 3 normal{q}(end) normal{q}(end)+2 normal{q}(end)+3];
    Nblock = length(blk);
    for k = 1:Nblock
        subplot(Ngroup,Nblock,(q-1)*Nblock + k); hold on
        plot([0 1],[0 0],'k')
        plot([0 0],[0 1],'k')
        for i = 1:Nfreq
            plot([0 mean(thetaOpt_mu{q}(1,blk(k),i,:),4)],[0 mean(thetaOpt_mu{q}(2,blk(k),i,:),4)],'LineWidth',1.5,'Color',map1(i,:))
            plot([0 mean(thetaOpt_mu{q}(3,blk(k),i,:),4)],[0 mean(thetaOpt_mu{q}(4,blk(k),i,:),4)],'LineWidth',1.5,'Color',map2(i,:))
        end
        axis([-0.45 1.2 -0.45 1.2])
        axis square
        xticks([])
        yticks([])
        if k == 1
            ylabel(groupNames{q})
        end
        if q == 1
            title(labels{k})
        end
    end
end

print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/vectors','-dpdf','-painters')

%% plot matrices as lines (normal blocks)
col = copper;
col = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

figure(2); clf
for q = 1:Ngroup
    mu = permute(thetaOpt_mu{q},[4 3 2 1]);
    se = permute(thetaOpt_se{q},[4 3 2 1]);
    Nblock = length(normal{q});
    totalTrials = Nblock*Ntrial;
    ticks = 1:5:totalTrials;
    ticks = ticks([1 2 end]);
    
    subplot(2, 3, q); hold on
    plot([1 totalTrials], [0 0], 'k')
    for k = 1:Nblock
%         if k <= Nblock
            block = normal{q}(k);
%         elseif k == Nblock + 1
%             block = flip{q}(1);
%         elseif k == Nblock + 2
%             block = flip{q}(2);
%         end
        plotIdx = Ntrial*(k-1)+(1:5);
        if k > 2
            plot([plotIdx(1)-0.5 plotIdx(1)-0.5],[-0.2 1],'Color',[0.8 0.8 0.8])
        end
        for i = 1:Nfreq
            s = shadedErrorBar(plotIdx, mu(:,i,block,1), se(:,i,block,1));
            editErrorBar(s,col(i,:),1);
        end
    end
    title(groupNames{q})
    axis([1 totalTrials -0.2 1])
    xticks(ticks)
    xticklabels({'Baseline','Early','Late'})
    set(gca,'TickDir','out')
    if q == 1
        ylabel('X --> X')
    end
    
    subplot(2, 3, q+3); hold on
    plot([1 totalTrials], [0 0], 'k')
    for k = 1:Nblock
%         if k <= Nblock
            block = normal{q}(k);
%         elseif k == Nblock + 1
%             block = flip{q}(1);
%         elseif k == Nblock + 2
%             block = flip{q}(2);
%         end
        
        plotIdx = Ntrial*(k-1)+(1:5);
        if k > 2
            plot([plotIdx(1)-0.5 plotIdx(1)-0.5],[-0.2 1],'Color',[0.8 0.8 0.8])
        end
        for i = 1:Nfreq
            s = shadedErrorBar(plotIdx, mu(:,i,block,4), se(:,i,block,4));
            editErrorBar(s,col(i,:),1);
        end
    end
    axis([1 totalTrials -0.2 1])
    xticks(ticks)
    xticklabels({'Baseline','Early','Late'})
    set(gca,'TickDir','out')
    if q == 1
        ylabel('Y --> Y')
    end
end

%%
labels = {'Late','Flip 1','Flip 2'};
figure(3); clf; 
for q = 1:Ngroup
    mu = permute(thetaOpt_mu{q},[4 3 2 1]);
    se = permute(thetaOpt_se{q},[4 3 2 1]);
    Nblock = length(flip{q})+1;
    totalTrials = Nblock*Ntrial;
    ticks = 1:5:totalTrials;
    
    subplot(2, 3, q); hold on
    plot([1 totalTrials], [0 0], 'k')
    for k = 1:Nblock
        if k == 1
            block = normal{q}(end);
        else
            block = flip{q}(k-1);
        end
        
        plotIdx = Ntrial*(k-1)+(1:5);
        for i = 1:Nfreq
            s = shadedErrorBar(plotIdx, mu(:,i,block,1), se(:,i,block,1));
            editErrorBar(s,col(i,:),1);
        end
    end
    title(groupNames{q})
    axis([1 totalTrials -0.5 1])
    xticks(ticks)
    xticklabels({'Baseline','Early','Late'})
    set(gca,'TickDir','out')
    if q == 1
        ylabel('X --> X')
    end
    
    subplot(2, 3, q+3); hold on
    plot([1 totalTrials], [0 0], 'k')
    for k = 1:Nblock
        if k == 1
            block = normal{q}(end);
        else
            block = flip{q}(k-1);
        end
        
        plotIdx = Ntrial*(k-1)+(1:5);
        for i = 1:Nfreq
            s = shadedErrorBar(plotIdx, mu(:,i,block,4), se(:,i,block,4));
            editErrorBar(s,col(i,:),1);
        end
    end
    axis([1 totalTrials -0.5 1])
    xticks(ticks)
    xticklabels({'Baseline','Early','Late'})
    set(gca,'TickDir','out')
    if q == 1
        ylabel('Y --> Y')
    end
end

%%
labels = {'Late','Flip 1','Flip 2'};
figure(4); clf; 
for q = 1:Ngroup
%     mu = permute(thetaOpt_mu{q},[4 3 2 1]);
%     se = permute(thetaOpt_se{q},[4 3 2 1]);
    Nblock = length(flip{q});
    totalTrials = Nblock*Ntrial;
    ticks = 1:5:totalTrials;
    
    percent = mean(abs(thetaOpt{q}(:,end-2,:,:,:)./thetaOpt{q}(:,end-4,:,:,:)),4);
    Nsubj = size(percent,5);
    mu = permute(mean(percent,5,'omitnan'),[3 1 2 4]);
    se = permute(std(percent,[],5,'omitnan'),[3 1 2 4])./sqrt(Nsubj);
    
    subplot(2, 3, q); hold on
    plot([1 totalTrials], [0 0], 'k')
%     for k = 1:Nblock
%         block = flip{q}(k);
        
%         plotIdx = Ntrial*(k-1)+(1:5);
        for i = 1:Nfreq
            errorbar(1,mu(i,1),se(i,1),'.','MarkerSize',20,'Color',col(i,:))
%             s = shadedErrorBar(plotIdx, mu(:,1), se(:,1));
%             editErrorBar(s,col(i,:),1);
        end
%     end
    title(groupNames{q})
    axis([0.5 1.5 -2 2])
    xticks(ticks)
    xticklabels({'Baseline','Early','Late'})
    set(gca,'TickDir','out')
    if q == 1
        ylabel('X --> X')
    end
    
    subplot(2, 3, q+3); hold on
    plot([1 totalTrials], [0 0], 'k')
%     for k = 1:Nblock
%         block = flip{q}(k);
        
%         plotIdx = Ntrial*(k-1)+(1:5);
        for i = 1:Nfreq
            errorbar(1,mu(i,4),se(i,4),'.','MarkerSize',20,'Color',col(i,:))
%             s = shadedErrorBar(plotIdx, mu(:,4), se(:,4));
%             editErrorBar(s,col(i,:),1);
        end
%     end
    axis([0.5 1.5 -2 2])
    xticks(ticks)
    xticklabels({'Baseline','Early','Late'})
    set(gca,'TickDir','out')
    if q == 1
        ylabel('Y --> Y')
    end
end
%% plot matrices as lines (dark blocks)
col = copper;
col = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

figure(4); clf
for q = 1:Ngroup
    mu = permute(thetaOpt_mu{q},[4 3 2 1]);
    se = permute(thetaOpt_se{q},[4 3 2 1]);
    Nblock = length(dark{q});
    totalTrials = Nblock*Ntrial;
    ticks = 1:5:totalTrials;
    ticks = ticks([1 2 end-1 end]);
    
    subplot(2, 3, q); hold on
    plot([1 totalTrials], [0 0], 'k')
    for k = 1:Nblock
        for i = 1:Nfreq
            s = shadedErrorBar(Ntrial*(k-1)+(1:5), mu(:,i,dark{q}(k),1), se(:,i,dark{q}(k),1));
            editErrorBar(s,col(i,:),1);
        end
    end
    title(groupNames{q})
    axis([1 totalTrials -0.5 1])
    xticks(ticks)
    xticklabels({'Baseline','Early','Late','Flip'})
    if q == 1
        ylabel('X --> X')
    end
    
    subplot(2, 3, q+3); hold on
    plot([1 totalTrials], [0 0], 'k')
    Nblock = length(dark{q});
    for k = 1:Nblock
        for i = 1:Nfreq
            s = shadedErrorBar(Ntrial*(k-1)+(1:5), mu(:,i,dark{q}(k),4), se(:,i,dark{q}(k),4));
            editErrorBar(s,col(i,:),1);
        end
    end
    axis([1 totalTrials -0.5 1])
    xticks(ticks)
    xticklabels({'Baseline','Early','Late','Flip'})
    if q == 1
        ylabel('Y --> Y')
    end
end

%% plot habit
col = [180 180 0
        0 191 255
        255 99 71]./255;
names = {'Early','Flip 1','Flip 2','Flip (dark)'};
    
figure(5); clf
for k = 1:4
    subplot(1,4,k); hold on
    plot([0.5 6.5], [0 0], 'k', 'HandleVisibility', 'off')
    for q = 1:Ngroup
        if k == 1
            block = 3;
        elseif k == 2
            block = size(thetaOpt{q},2)-2;
        elseif k == 3
            block = size(thetaOpt{q},2)-1;
        elseif k == 4
            block = size(thetaOpt{q},2);
        end
        
        h = squeeze(mean(thetaOpt{q}(1,block,:,:,:),4));
        mu = mean(h,2);
        se = std(h,[],2)./sqrt(size(h,2));
        
        s = shadedErrorBar(1:Nfreq,mu,se);
        editErrorBar(s,col(q,:),1);
    end
    axis([0.5 6.5 -0.5 0.4])
    xticks(1:6)
    xlabel('Frequency')
    ylabel('Gain')
    title(names{k})
    set(gca, 'TickDir', 'out')
end
legend(groupNames)

%% plot habit
% figure(5); clf
% for q = 1:Ngroup
%     n = length(normal{q});
%     f = length(flip{q});
%     for k = 1:Nfreq
%         subplot(3,Nfreq,(q-1)*6+k); hold on
%         plot([0 Nblock],[0 0],'k','LineWidth',1)
%         plot(1:n,squeeze(thetaOpt{q}(response,normal{q},k,:)),'Color',[0 0 0 0.6])
%         plot(n+1:n+f,squeeze(thetaOpt{q}(response,flip{q},k,:)),'.','Color',[0 0 0 0.6],'MarkerSize',15)
%         plot(1:n,mat{q}(response,normal{q},k)','Color',col(2,:),'LineWidth',3)
%         plot(n+1:n+f,mat{q}(response,flip{q},k)','.','Color',col(2,:),'MarkerSize',30)
%         axis([1 maxNormal+1 -0.6 1])
%         xticks([])
%         yticks([])
%         if k == 1
%             title('Low freq')
%             xlabel('Blocks')
%             ylabel('Gain (response axis)')
%             yticks(-0.6:0.2:1)
%         elseif k == Nfreq
%             title('High freq')
%         end
%     end
% end
% 
% figure(6); clf
% for q = 1:Ngroup
%     subplot(2,3,q); hold on
%     plot([1 Nfreq],[0 0],'k','LineWidth',1)
%     plot(squeeze(thetaOpt{q}(response,normal{q}(end),:,:)),'Color',[0 0 0 0.6])
%     plot(squeeze(mat{q}(response,normal{q}(end),:)),'Color',col(2,:),'LineWidth',3)
%     xticks([1 6])
%     xticklabels({'Low Freq','High Freq'})
%     yticks(-0.9:0.3:0.9)
%     if q == 1
%         ylabel('Gain (before flip; response axis)')
%     end
%     axis([1 6 -0.7 0.9])
%     title('Before flip')
%     
%     subplot(2,3,q+3); hold on
%     plot([1 Nfreq],[0 0],'k','LineWidth',1)
%     plot(squeeze(thetaOpt{q}(response,flip{q},:,:)),'Color',[0 0 0 0.6])
%     plot(squeeze(mat{q}(response,flip{q},:)),'Color',col(2,:),'LineWidth',3)
%     xticks([1 6])
%     xticklabels({'Low Freq','High Freq'})
%     yticks(-0.9:0.3:0.9)
%     if q == 1
%         ylabel('Gain (after flip; response axis)')
%     end
%     axis([1 6 -0.7 0.9])
% end

%% plot strength of habit within subject
figure(7); clf
for q = 1:Ngroup
    for i = 1:Nfreq
        habit = squeeze([thetaOpt{q}(response,normal{q}(end),i,:) thetaOpt{q}(response,flip{q},i,:)]);
        subplot(Ngroup,Nfreq,(q-1)*Nfreq+i); hold on
        plot([1 2],[0 0],'k')
        plot(habit,'Color',[0 0 0 0.6])
        plot(mean(habit,2),'.k','MarkerSize',30)
        xticks([])
        ylim([-0.5 1])
        if q == 1
            if i == 1
                ylabel('Gain (2 day)')
                title('Low frequency')
            elseif i == 6
                title('High frequency')
            end
        elseif q == 2 && i == 1
            ylabel('Gain (5 day)')
        elseif q == 3
            xticks([1 2])
            xticklabels({'Before flip','After flip'})
            if i == 1
                ylabel('Gain (10 day)')
            end
        end
    end
end

%% compare habit across groups
figure(8); clf
for i = 1:Nfreq
    subplot(1,6,i); hold on
    plot([1 3],[0 0],'k')
    for q = 1:Ngroup
        plot(q,squeeze(thetaOpt{q}(response,flip{q},i,:)),'.k','MarkerSize',15)
        plot(q,mean(thetaOpt{q}(response,flip{q},i,:),4),'.','Color',col(2,:),'MarkerSize',30)
    end
    xticks(1:3)
    xticklabels({'2 day','5 day','10 day'})
    ylim([-0.6 0.6])
end

%% compare habit to baseline
figure(9); clf
for q = 1:Ngroup
    for i = 1:Nfreq
        habit = squeeze([thetaOpt{q}(response,normal{q}(1),i,:) thetaOpt{q}(response,flip{q},i,:)]);
        subplot(Ngroup,Nfreq,(q-1)*Nfreq+i); hold on
        plot([1 2],[0 0],'k')
        plot(habit,'Color',[0 0 0 0.6])
        plot(mean(habit,2),'.','Color',col(2,:),'MarkerSize',30)
        ylim([-.8 .8])
    end
end

%% plot proportion of movement in workspace vs null space
workTotal = squeeze(sum(abs(thetaOpt(workspace,:,:,:)),1));
nullTotal = squeeze(sum(abs(thetaOpt(nullspace,:,:,:)),1));
total = workTotal+nullTotal;
workProp = workTotal./total;

baseline_blocks = 1:2;
workBase = squeeze(sum(abs(thetaOpt([1 4],baseline_blocks,:,:)),1));
nullBase = squeeze(sum(abs(thetaOpt([2 3],baseline_blocks,:,:)),1));
totalBase = workBase + nullBase;
baseProp = workBase./totalBase;

workProp(baseline_blocks,:,:) = baseProp;
workPropMu = mean(workProp,3);

figure(11); clf
for k = 1:Nfreq
    subplot(2,Nfreq,k); hold on
    plot(workPropMu(normal,k),'LineWidth',2)
    for i = 1:Nsubj
        plot(workProp(normal,k,i),'Color',[col(1,:) 0.6])
    end
    if k == 1
        title('Normal blocks')
        ylabel('Proportion of movements in workspace')
    end
    axis([1 length(normal) 0 1])

    subplot(2,Nfreq,Nfreq+k); hold on
    plot(workPropMu(special,k),'LineWidth',2)
    for i = 1:Nsubj
        plot(workProp(special,k,i),'Color',[col(1,:) 0.6])
    end
    if k == 1
        title('Special blocks')
        xlabel('Low Freq')
        ylabel('Proportion of movements in workspace')
    elseif k == Nfreq
        xlabel('High freq')
    end
    axis([1 length(normal) 0 1])
end

%% gain matrices
figure(3); clf
for i = 1:Nfreq
    subplot(1,Nfreq,i)
    imagesc(mat(:,index,i),clims)
    colormap(map)
    set(gca,'Xtick',[],'Ytick',[])
    if i == 1
        title('Low freq')
        yticks(1:4)
        yticklabels(labels)
    elseif i == Nfreq
        title('High freq')
        if remove == 0
            xticks(special)
            xticklabels(graph_name(special))
        end
    end
end
%% plot gain matrices as 2x2
subj = 6;
rMat = rotMat(:,:,:,:,subj);

% rMat = mean(rotMat,5);

figure(4); clf
for i = 1:Nfreq
    for k = 1:length(index)
        subplot(Nfreq,length(index),(i-1)*length(index)+k)
        imagesc(rMat(:,:,i,k),clims)
        colormap(map)
        pbaspect([1 1 1])
        set(gca,'Xtick',[],'Ytick',[])
        if i == 1
            title(graph_name(k))
        end
        if k == 1 && i == 1
            ylabel('Low freq')
        elseif k == 1 && i == Nfreq
            ylabel('High freq')
        end
    end
end

%% display results of fitting process
g = 1;
subj = 1;
blockIdx = 1;
freq = 1;

% slope is the line which all data for a given frequency projects onto
slope = template(:,freq,subj,g);
for i = 1:numel(slope)
    slope(i) = slope(i)*(2/norm(slope(i)));
end

% reconstruct fitted phasor
template2 = template(:,:,subj,g);
template2 = [repmat(template2(1,:),[2 1]); repmat(template2(2,:),[2 1])];
phasor = squeeze(thetaOpt(:,blockIdx,:,subj,g)).*template2;
dat = squeeze(cplx(freq,:,:,subj,g)); % extract data to plot

gblocks = [1 3 2 4];
plot_subj(data{subj}.(block_name{blockIdx}).phasors.(output))
for i = 1:4
    subplot(2,2,gblocks(i)); hold on
    plot(phasor(i,:),'-ok','LineWidth',1.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none')
    plot([phasor(i,freq) cplx(freq,i,blockIdx,subj,g)],'b','LineWidth',2)
    if i <=2
        plot(dat([1 2],:),'.b','MarkerSize',10)
        if thetaOpt(i,blockIdx,freq) >= 0
            plot([0 real(slope(1))],[0 imag(slope(1))],'b','LineWidth',2)
        else
            plot([0 -real(slope(1))],[0 -imag(slope(1))],'b','LineWidth',2)
        end
    else
        plot(dat([3 4],:),'.b','MarkerSize',10)
        if thetaOpt(i,blockIdx,freq) >= 0
            plot([0 real(slope(2))],[0 imag(slope(2))],'b','LineWidth',2)
        else
            plot([0 -real(slope(2))],[0 -imag(slope(2))],'b','LineWidth',2)
        end
    end
    axis([-1.25 1.25 -1.25 1.25])
    pbaspect([1 1 1])
end

%% graph average fitted phasors across subjects
col1 = [1 0.9 0.3];
col2 = [1 0 0];
col = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2),Nfreq)', linspace(col1(3),col2(3),Nfreq)'];
gblocks = [1 3 2 4];

group = 1;
blockIdx = 6;

phasorMat = repmat(permute(phasor(:,:,group),[1 3 2]),[1 4]);
phasorLines = permute(thetaOpt(:,blockIdx,:,:,group),[3 1 4 2]).*phasorMat;
phasorLines = mean(phasorLines,3);

dat_all = permute(squeeze(cplx(:,:,blockIdx,:,group)),[3 1 2]);
dat = squeeze(mean(dat_all,1));

figure(6); clf
for k = 1:4
    subplot(2,2,gblocks(k)); hold on
    for i = 1:Nfreq
        plot(dat_all(:,i,k),'.','Color',col(i,:),'MarkerSize',10)
        if i > 1
            plot(dat(i-1:i,k),'-o','Color',col(i,:),'LineWidth',1,'MarkerFaceColor',col(i-1,:),'MarkerEdgeColor','none')
        end
    end
    plot(phasorLines(:,k),'-ok','LineWidth',2,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none')
    plot([-1.5 1.5],[0 0],'k','LineWidth',1)
    plot([0 0],[-1.5 1.5],'k','LineWidth',1)
    axis([-1.25 1.25 -1.25 1.25])
    pbaspect([1 1 1])
end

end