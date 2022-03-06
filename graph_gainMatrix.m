% graph_gainMatrix plots the direction and amplitude of cursor movements in
% response to target movement at different frequencies
% 
%   data: structure containing all data
%   block_name: name of blocks
%   blockType: identifies whether the block has visual feedback (=1), no 
%       visual feedback (=2), or the flipped mapping (=3)

function graph_gainMatrix(data, block_name, blockType)

output = 'cursor';
groups = {'day2','day5','day10'};
groupNames = {'2-day','5-day','10-day'};
names = {'x_x','x_y','y_x','y_y'};

Ngroup = length(groups);

for q = 1:Ngroup % loop over groups
    Nsubj = length(data.(groups{q}));
    blk = block_name.(groups{q});
    Nblock = length(blk);
    Nfreq = length(data.(groups{q}){1}.(blk{1}).freqX);
    Ntrial = 5;
    
    normal{q} = find(blockType.(groups{q}) == 1);
    dark{q} = find(blockType.(groups{q}) == 2);
    
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
            for j = 1:2
            
                if j == 1 % average baseline phasor
                    comparison = mean(cplx{q}(i,:,[1 4],1,p),2,'omitnan');
                    blocks = 1:normal{q}(end);
                    idx = 1:normal{q}(end)*2*Ntrial;
                else % average late learning phasor
                    comparison = mean(cplx{q}(i,:,[1 4],normal{q}(end),p),2,'omitnan');
                    blocks = (normal{q}(end)+1):(normal{q}(end)+4);
                    idx = (normal{q}(end)*2*Ntrial + 1):((normal{q}(end)+4)*2*Ntrial);
                end
                
                % x-component of xx and yy phasor
                x = cos(angle(comparison));
                
                % y-component of xx and yy phasor
                y = sin(angle(comparison));
                
                % project data onto baseline x phasor
                dat = cplx{q}(i,:,[1 2],blocks,p); % extract data for x-target to x-/y-hand
                phasor = reshape(dat, [1 numel(dat)]);
                num = [real(phasor); imag(phasor)]; % separate real and imaginary parts of phasors
                unitVec = [x(1); y(1)]; % unit vector with phase equal to baseline vector
                gainX{q}(idx,i,p) = dot(num, repmat(unitVec, [1 numel(dat)])); % project num onto unitVec
                
                % project data onto baseline y phasor
                dat = cplx{q}(i,:,[3 4],blocks,p); % extract data for x-target to x-/y-hand
                phasor = reshape(dat, [1 numel(dat)]);
                num = [real(phasor); imag(phasor)]; % separate real and imaginary parts of phasors
                unitVec = [x(2); y(2)]; % unit vector with phase equal to baseline vector
                gainY{q}(idx,i,p) = dot(num, repmat(unitVec, [1 numel(dat)])); % project num onto unitVec
                
            end
        end
    end
    
    % combine gainX and gainY
    thetaOpt{q} = [reshape(gainX{q},[Ntrial 2 Nblock Nfreq Nsubj]) ...
        reshape(gainY{q},[Ntrial 2 Nblock Nfreq Nsubj])];
    thetaOpt{q} = permute(thetaOpt{q}, [2 3 4 1 5]);
    
    % average and standard error
    thetaOpt_mu{q} = mean(thetaOpt{q},5,'omitnan');
    thetaOpt_se{q} = std(thetaOpt{q},[],5,'omitnan')./sqrt(Nsubj);
end

% store data for analysis in R
blockNum = [3 5; 5 11; 11 21];
y = [];
for k = 1:Ngroup
    for i = 1:2        
        t = mean(thetaOpt{k}(1,blockNum(k,i),:,:,:),4,'omitnan');
        y = [y; t(:)];
    end
end
g(1:156,1) = "2-day";
g(157:324,1) = "5-day";
g(325:384,1) = "10-day";
b([1:78 157:240 325:354],1) = "before";
b([79:156 241:324 355:384],1) = "after";
frequency = repmat((1:Nfreq)',[64 1]);
s1 = repmat(1:13,[Nfreq 2]);
s2 = repmat(14:27,[Nfreq 2]);
s3 = repmat(28:32,[Nfreq 2]);
subject = [s1(:); s2(:); s3(:)];
T = table(g, b, frequency, subject, y, 'VariableNames', {'group','block','frequency','subject','gain'});
% writetable(T,'C:/Users/Chris/Documents/R/habit/data/gain_skill.csv')

% for plotting vectors
col1 = [0 128 0]/255;
col2 = [152 251 152]/255;
map1 = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2),Nfreq)', linspace(col1(3),col2(3),Nfreq)'];

col1 = [128 0 128]/255;
col2 = [230 230 250]/255;
map2 = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2),Nfreq)', linspace(col1(3),col2(3),Nfreq)'];

%% plot vectors

labels = {'Baseline','Early','Late',"Flip 1","Flip 2"};
figure(4); clf
for q = 1:Ngroup
    blk = [1 3 normal{q}(end) normal{q}(end)+2 normal{q}(end)+3];
    Nblock = length(blk);
    for k = 1:Nblock
        subplot(Ngroup,Nblock,(q-1)*Nblock + k); hold on
        plot([0 0.5],[0 0],'k')
        plot([0 0],[0 0.5],'k')
        for i = 1:Nfreq
            plot([0 mean(thetaOpt_mu{q}(1,blk(k),i,:),4)],[0 mean(thetaOpt_mu{q}(2,blk(k),i,:),4)],'LineWidth',1.5,'Color',map1(i,:))
            plot([0 mean(thetaOpt_mu{q}(3,blk(k),i,:),4)],[0 mean(thetaOpt_mu{q}(4,blk(k),i,:),4)],'LineWidth',1.5,'Color',map2(i,:))
        end
        axis([-0.35 0.8 -0.35 0.8])
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

% save figure for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/vectors','-dpdf','-painters')

%% plot matrices as lines (normal blocks)
offset = [0 20 55];

f = figure(5); clf
set(f,'Position',[200 200 380 250]);

for j = 1:2
    subplot(2,1,j); hold on
    for q = 1:Ngroup
        mu = permute(thetaOpt_mu{q},[4 3 2 1]);
        se = permute(thetaOpt_se{q},[4 3 2 1]);
        Nblock = length(normal{q});
        totalTrials = Nblock*Ntrial;
        
        plot([offset(q)+1 offset(q)+totalTrials], [0 0], 'k')
        plot([offset(q)+1 offset(q)+totalTrials], [-.2 -.2], 'k')
        if q > 1
            plot([offset(q)+1 offset(q)+1], [-.2 1], 'k')
        end
        for k = 1:Nblock
            block = normal{q}(k);
            plotIdx = Ntrial*(k-1)+(1:5) + offset(q);
            if k > 2
                plot([plotIdx(1)-0.5 plotIdx(1)-0.5],[-0.2 1],'Color',[0.8 0.8 0.8])
            end
            for i = 1:Nfreq
                if j == 1
                    s = shadedErrorBar(plotIdx, mu(:,i,block,1), se(:,i,block,1));
                else
                    s = shadedErrorBar(plotIdx, mu(:,i,block,4), se(:,i,block,4));
                end
                editErrorBar(s,map1(i,:),1);
            end
        end
    end
    axis([1 offset(3)+totalTrials -0.2 1])
    yticks(0:0.25:1)
    set(gca,'TickDir','out','Xcolor','none')
    if j == 1
        ylabel('X --> X')
    else
        ylabel('Y --> Y')
    end
end

% save figure for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/gains','-dpdf','-painters')

%%

f = figure(6); clf
set(f,'Position',[200 200 380 250]);

for j = 1:2
    subplot(2,1,j); hold on
    for q = 1:Ngroup
        mu = permute(thetaOpt_mu{q},[4 3 2 1]);
        se = permute(thetaOpt_se{q},[4 3 2 1]);
        Nblock = length(dark{q})-1;
        totalTrials = Nblock*Ntrial;
        
        plot([offset(q)+1 offset(q)+totalTrials], [0 0], 'k')
        plot([offset(q)+1 offset(q)+totalTrials], [-.3 -.3], 'k')
        if q > 1
            plot([offset(q)+1 offset(q)+1], [-.4 1.25], 'k')
        end
        for k = 1:Nblock
            block = dark{q}(k);
            plotIdx = Ntrial*(k-1)+(1:5) + offset(q);
            if k > 2
                plot([plotIdx(1)-0.5 plotIdx(1)-0.5],[-0.4 1.25],'Color',[0.8 0.8 0.8])
            end
            for i = 1:Nfreq
                if j == 1
                    s = shadedErrorBar(plotIdx, mu(:,i,block,1), se(:,i,block,1));
                else
                    s = shadedErrorBar(plotIdx, mu(:,i,block,4), se(:,i,block,4));
                end
                editErrorBar(s,map1(i,:),1);
            end
        end
    end
    axis([1 offset(3)+totalTrials -0.3 1.25])
    yticks(-0.5:0.5:1)
    set(gca,'TickDir','out','Xcolor','none')
    if j == 1
        ylabel('X --> X')
    else
        ylabel('Y --> Y')
    end
end

% save figure for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/gains_dark','-dpdf','-painters')

end