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
    
    dark{q} = find(contains(graph_name.(groups{q}),'(D)'));
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
    
    % convert gain matrices from cursor to hand configuration
%     thetaOpt{q}(:,3:end,:,:,:) = thetaOpt{q}([2 1 4 3],3:end,:,:,:);
    
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

% save flipSign flipSign

y = [];
for k = 1:Ngroup
    for i = 1:2
        if i == 1
            t = mean(thetaOpt{k}(4,1,:,:,:),4,'omitnan');
        else
            t = mean(thetaOpt{k}(1,end-4,:,:,:),4,'omitnan');
        end
        y = [y; t(:)];
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
            if k == 1
                plot([0 mean(thetaOpt_mu{q}(1,blk(k),i,:),4)],[0 mean(thetaOpt_mu{q}(2,blk(k),i,:),4)],'LineWidth',1.5,'Color',map1(i,:))
                plot([0 mean(thetaOpt_mu{q}(3,blk(k),i,:),4)],[0 mean(thetaOpt_mu{q}(4,blk(k),i,:),4)],'LineWidth',1.5,'Color',map2(i,:))
            else
                plot([0 mean(thetaOpt_mu{q}(2,blk(k),i,:),4)],[0 mean(thetaOpt_mu{q}(1,blk(k),i,:),4)],'LineWidth',1.5,'Color',map1(i,:))
                plot([0 mean(thetaOpt_mu{q}(4,blk(k),i,:),4)],[0 mean(thetaOpt_mu{q}(3,blk(k),i,:),4)],'LineWidth',1.5,'Color',map2(i,:))
            end
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

% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/vectors','-dpdf','-painters')

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
        block = normal{q}(k);
        plotIdx = Ntrial*(k-1)+(1:5);
        if k > 2
            plot([plotIdx(1)-0.5 plotIdx(1)-0.5],[-0.2 1],'Color',[0.8 0.8 0.8])
        end
        for i = 1:Nfreq
            if k == 1
                s = shadedErrorBar(plotIdx, mu(:,i,block,1), se(:,i,block,1));
                editErrorBar(s,col(i,:),1);
            else
                s = shadedErrorBar(plotIdx, mu(:,i,block,4), se(:,i,block,4));
                editErrorBar(s,col(i,:),1);
            end
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
        block = normal{q}(k);
        
        plotIdx = Ntrial*(k-1)+(1:5);
        if k > 2
            plot([plotIdx(1)-0.5 plotIdx(1)-0.5],[-0.2 1],'Color',[0.8 0.8 0.8])
        end
        for i = 1:Nfreq
            if k == 1
                s = shadedErrorBar(plotIdx, mu(:,i,block,4), se(:,i,block,4));
                editErrorBar(s,col(i,:),1);
            else
                s = shadedErrorBar(plotIdx, mu(:,i,block,1), se(:,i,block,1));
                editErrorBar(s,col(i,:),1);
            end
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
            if k == 1
                s = shadedErrorBar(plotIdx, mu(:,i,block,1), se(:,i,block,1));
                editErrorBar(s,col(i,:),1);
            else
                s = shadedErrorBar(plotIdx, mu(:,i,block,4), se(:,i,block,4));
                editErrorBar(s,col(i,:),1);
            end
        end
    end
    title(groupNames{q})
    axis([1 totalTrials -0.5 1])
    xticks(ticks)
    xticklabels(labels)
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
            if k == 1
                s = shadedErrorBar(plotIdx, mu(:,i,block,4), se(:,i,block,4));
                editErrorBar(s,col(i,:),1);
            else
                s = shadedErrorBar(plotIdx, mu(:,i,block,1), se(:,i,block,1));
                editErrorBar(s,col(i,:),1);
            end
        end
    end
    axis([1 totalTrials -0.5 1])
    xticks(ticks)
    xticklabels(labels)
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
    ticks = ticks([1 2 end]);
    
    subplot(2, 3, q); hold on
    plot([1 totalTrials], [0 0], 'k')
    for k = 1:Nblock
        for i = 1:Nfreq
            if k == 1
                s = shadedErrorBar(Ntrial*(k-1)+(1:5), mu(:,i,dark{q}(k),1), se(:,i,dark{q}(k),1));
                editErrorBar(s,col(i,:),1);
            else
                s = shadedErrorBar(Ntrial*(k-1)+(1:5), mu(:,i,dark{q}(k),4), se(:,i,dark{q}(k),4));
                editErrorBar(s,col(i,:),1);
            end
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
            if k == 1
                s = shadedErrorBar(Ntrial*(k-1)+(1:5), mu(:,i,dark{q}(k),4), se(:,i,dark{q}(k),4));
                editErrorBar(s,col(i,:),1);
            else
                s = shadedErrorBar(Ntrial*(k-1)+(1:5), mu(:,i,dark{q}(k),1), se(:,i,dark{q}(k),1));
                editErrorBar(s,col(i,:),1);
            end
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

end