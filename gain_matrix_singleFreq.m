output = 'Lhand';
names = {'x_x_all','x_y_all','y_x_all','y_y_all'};
Nblock = length(block_name);
Nsubj = length(data)-1;
Nfreq = length(data{end}.freqX);
Ngroups = 1;
paramsInit = zeros([2 Nblock]);
special = find(contains(graph_name,'('));

% set remove = 1 to remove special blocks
index = 1:Nblock;
normal = index;
normal(special) = [];
remove = 0;
if remove
    index = normal;
end

if strcmp(output,'cursor')
    labels = {'X_{T}X_{O} (response)','X_{T}Y_{O}','Y_{T}X_{O}','Y_{T}Y_{O} (response)'};
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

for p = 1:Nsubj % loop over subjects
    for k = 1:Nblock
        for i = 1:4
            cplx(:,i,k,p) = mean(data{p}.(block_name{k}).phasors.(output).(names{i}).ratio,2); % put all complex ratios in cplx
        end
    end
    
    for i = 1:Nfreq
        x = cos(angle(cplx(i,[1 4],1,p)));
        y = sin(angle(cplx(i,[1 4],1,p)));
        template(:,i,p) = x+y*1j; % generate template of gain = 1 and phase same as baseline
        error = @(params) scale(params,cplx(i,[1 2],:,p),template(1,i,p));
        opt1(:,:,i,p) = fmincon(error,paramsInit); % fit gains for x frequencies
        error = @(params) scale(params,cplx(i,[3 4],:,p),template(2,i,p));
        opt2(:,:,i,p) = fmincon(error,paramsInit); % fit gains for y frequencies
    end
end

% combine opt1 and opt2
thetaOpt = [reshape(opt1,[2 Nblock Nfreq Nsubj Ngroups]); reshape(opt2,[2 Nblock Nfreq Nsubj Ngroups])];

% shape thetaOpt into gain matrix format
rotMat = reshape(thetaOpt,[2 2 Nblock Nfreq Nsubj Ngroups]);
rotMat = permute(rotMat,[1 2 4 3 5 6]);

% for averaging across subjects
mat = squeeze(mean(thetaOpt,4));

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

%% plot vectors and gain matrices
gblocks = 1:3;
figure(2); clf
for k = 1:length(gblocks)
    subplot(1,length(gblocks),k); hold on
    plot([0 1],[0 0],'k')
    plot([0 0],[0 1],'k')
    for i = 1:Nfreq
        plot([0 mat(1,gblocks(k),i)],[0 mat(2,gblocks(k),i)],'LineWidth',1.5,'Color',map1(i,:))
        plot([0 mat(3,gblocks(k),i)],[0 mat(4,gblocks(k),i)],'LineWidth',1.5,'Color',map2(i,:))
        axis([-0.45 1.2 -0.45 1.2])
        axis square
    end
    title(graph_name{gblocks(k)})
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gain matrices
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

%% plot matrices as lines
ymin = min(min(min(mat)));
ymax = max(max(max(mat)));

figure(5); clf
for k = 1:Nfreq
    subplot(2,Nfreq,k); hold on
    plot(mat(workspace,normal,k)','LineWidth',2)
    for i = 1:Nsubj
        for j = 1:2
            plot(thetaOpt(workspace(j),normal,k,i),'Color',[col(j,:) 0.6])
        end
    end
    plot([0 Nblock+1],[0 0],'--k','LineWidth',1)
    axis([1 length(normal) ymin ymax])
    xticks([])
    if k == 1
        title('Normal blocks')
        ylabel('Gain (workspace)')
        xlabel('Low freq')
    elseif k == Nfreq
        xlabel('High freq')
        legend(labels(workspace))
    end
    
    subplot(2,Nfreq,k+Nfreq); hold on
    set(gca,'ColorOrderIndex',3)
    plot(mat(nullspace,normal,k)','LineWidth',2)
    for i = 1:Nsubj
        for j = 1:2
            plot(thetaOpt(nullspace(j),normal,k,i),'Color',[col(j+2,:) 0.6])
        end
    end
    plot([0 Nblock+1],[0 0],'--k','LineWidth',1)
    axis([1 length(normal) ymin ymax])
    xticks([])
    if k == 1
        ylabel('Gain (null space)')
        xlabel('Low freq')
    elseif k == Nfreq
        xlabel('High freq')
        legend(labels(nullspace))
    end
end

figure(6); clf
for k = 1:Nfreq
    subplot(2,Nfreq,k); hold on
    plot(mat(workspace,special,k)','LineWidth',2)
    for i = 1:Nsubj
        for j = 1:2
            plot(thetaOpt(workspace(j),special,k,i),'Color',[col(j,:) 0.6])
        end
    end
    plot([0 Nblock+1],[0 0],'--k','LineWidth',1)
    axis([1 length(special) ymin ymax])
    xticks([])
    if k == 1
        title('Special blocks')
        ylabel('Gain (workspace)')
        xlabel('Low freq')
    elseif k == Nfreq
        xlabel('High freq')
        legend(labels(workspace))
    end
    
    subplot(2,Nfreq,k+Nfreq); hold on
    set(gca,'ColorOrderIndex',3)
    plot(mat(nullspace,special,k)','LineWidth',2)
    for i = 1:Nsubj
        for j = 1:2
            plot(thetaOpt(nullspace(j),special,k,i),'Color',[col(j+2,:) 0.6])
        end
    end
    plot([0 Nblock+1],[0 0],'--k','LineWidth',1)
    axis([1 length(special) ymin ymax])
    xticks([])
    if k == 1
        ylabel('Gain (null space)')
        xlabel('Low freq')
    elseif k == Nfreq
        xlabel('High freq')
        legend(labels(nullspace))
    end
end

figure(7); clf
for k = 1:Nfreq
    subplot(2,Nfreq,k); hold on
    plot(mat(workspace,:,k)','LineWidth',2)
    for i = 1:Nsubj
        for j = 1:2
            plot(thetaOpt(workspace(j),:,k,i),'Color',[col(j,:) 0.6])
        end
    end
    plot([0 Nblock+1],[0 0],'--k','LineWidth',1)
    axis([1 size(mat,2) ymin ymax])
    xticks([])
    if k == 1
        title('All blocks')
        ylabel('Gain (workspace)')
        xlabel('Low freq')
    elseif k == Nfreq
        xlabel('High freq')
        legend(labels(workspace))
    end
    
    subplot(2,Nfreq,k+Nfreq); hold on
    set(gca,'ColorOrderIndex',3)
    plot(mat(nullspace,:,k)','LineWidth',2)
    for i = 1:Nsubj
        for j = 1:2
            plot(thetaOpt(nullspace(j),:,k,i),'Color',[col(j+2,:) 0.6])
        end
    end
    plot([0 Nblock+1],[0 0],'--k','LineWidth',1)
    axis([1 size(mat,2) ymin ymax])
    xticks([])
    if k == 1
        ylabel('Gain (null space)')
        xlabel('Low freq')
    elseif k == Nfreq
        xlabel('High freq')
        legend(labels(nullspace))
    end
end

%% plot habit
se = squeeze(std(thetaOpt(2,:,:,:),[],4))/sqrt(Nsubj);
n = length(normal);
s = length(special);
figure(8); clf
for k = 1:Nfreq
    subplot(1,Nfreq,k); hold on
    errorbar(1:n,mat(2,normal,k)',se(1:n,k),'-ko','LineWidth',2,'MarkerFaceColor','k')
    errorbar(n+1:n+s,mat(2,special,k)',se(n+1:end,k),'-ko','LineWidth',2,'MarkerFaceColor','k')
%     plot(1:length(normal),mat(2,normal,k)','-ko','LineWidth',2,'MarkerFaceColor','k')
%     plot(length(normal)+1:length(normal)+length(special),mat(2,special,k)','-ko','LineWidth',2,'MarkerFaceColor','k')
    for i = 1:Nsubj
        plot(1:length(normal),thetaOpt(2,normal,k,i),'Color',[0 0 0 0.6])
        plot(length(normal)+1:length(normal)+length(special),thetaOpt(2,special,k,i),'Color',[0 0 0 0.6])
    end
    plot([0 Nblock+1],[0 0],'--k','LineWidth',1)
    axis([1 size(mat,2) -0.7 1])
    xticks([])
    yticks(-0.6:0.3:0.9)
    if k == 1
        title('Low freq')
        xlabel('Blocks')
        ylabel('Gain (response axis)')
    elseif k == Nfreq
        title('High freq')
    end
end

figure(9); clf
subplot(1,2,1); hold on
plot([1 Nfreq],[0 0],'--k','LineWidth',1)
errorbar(squeeze(mat(2,6,:)),se(6,:),'-ko','LineWidth',2,'MarkerFaceColor','k')
% plot(squeeze(mat(2,5,:)),'-ko','LineWidth',2,'MarkerFaceColor','k')
for i = 1:Nsubj
    plot(squeeze(thetaOpt(2,6,:,i)),'Color',[0 0 0 0.6])
end
xticks([1 6])
xticklabels({'Low Freq','High Freq'})
yticks(-0.9:0.3:0.9)
ylabel('Gain (response axis)')
axis([1 6 -0.7 0.9])
title('Before flip')

subplot(1,2,2); hold on
plot([1 Nfreq],[0 0],'--k','LineWidth',1)
errorbar(squeeze(mat(2,7,:)),se(7,:),'-ko','LineWidth',2,'MarkerFaceColor','k')
% plot(squeeze(mat(2,6,:)),'-ko','LineWidth',2,'MarkerFaceColor','k')
for i = 1:Nsubj
    plot(squeeze(thetaOpt(2,7,:,i)),'Color',[0 0 0 0.6])
end
xticks([1 6])
xticklabels({'Low Freq','High Freq'})
yticks(-0.9:0.3:0.9)
axis([1 6 -0.7 0.9])
title('After flip')

%% plot effect of dual task
effect = thetaOpt(:,normal,:,:) - thetaOpt(:,special,:,:);
effectMu = mean(effect,4);

ymin = min(min(min(effectMu)));
ymax = max(max(max(effectMu)));

figure(9); clf
for k = 1:Nfreq
    subplot(2,Nfreq,k); hold on
    plot(effectMu(workspace,:,k)','LineWidth',2)
    for i = 1:Nsubj
        for j = 1:2
            plot(effect(workspace(j),:,k,i),'Color',[col(j,:) 0.6])
        end
    end
    plot([0 Nblock+1],[0 0],'--k','LineWidth',1)
    axis([1 length(normal) ymin ymax])
    xticks([])
    if k == 1
        title('Normal blocks')
        ylabel('Gain (workspace)')
        xlabel('Low freq')
    elseif k == Nfreq
        xlabel('High freq')
        legend(labels(workspace))
    end
    
    subplot(2,Nfreq,k+Nfreq); hold on
    set(gca,'ColorOrderIndex',3)
    plot(effectMu(nullspace,:,k)','LineWidth',2)
    for i = 1:Nsubj
        for j = 1:2
            plot(effect(nullspace(j),:,k,i),'Color',[col(j+2,:) 0.6])
        end
    end
    plot([0 Nblock+1],[0 0],'--k','LineWidth',1)
    axis([1 length(normal) ymin ymax])
    xticks([])
    if k == 1
        ylabel('Gain (null space)')
        xlabel('Low freq')
    elseif k == Nfreq
        xlabel('High freq')
        legend(labels(nullspace))
    end
end

%%
normalized = (thetaOpt([1 4 response],normal,:,:) - thetaOpt([1 4 response],special,:,:))./thetaOpt([1 4 response],normal,:,:);
baseline = mean(normalized(1:2,:,:,:),1);
effect_norm = squeeze([baseline(:,1,:,:) normalized(3,2:end,:,:)]);
effectMu_norm = mean(effect_norm,3);

figure(10); clf
plot(effectMu_norm)

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
colors = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2),Nfreq)', linspace(col1(3),col2(3),Nfreq)'];
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
        plot(dat_all(:,i,k),'.','Color',colors(i,:),'MarkerSize',10)
        if i > 1
            plot(dat(i-1:i,k),'-o','Color',colors(i,:),'LineWidth',1,'MarkerFaceColor',colors(i-1,:),'MarkerEdgeColor','none')
        end
    end
    plot(phasorLines(:,k),'-ok','LineWidth',2,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none')
    plot([-1.5 1.5],[0 0],'k','LineWidth',1)
    plot([0 0],[-1.5 1.5],'k','LineWidth',1)
    axis([-1.25 1.25 -1.25 1.25])
    pbaspect([1 1 1])
end

%%
function e = scale(params,cplx_data,template)
    phasors = params*template;
    e = (phasors - squeeze(cplx_data)).^2;
    for k = 1:numel(e)
        e(k) = norm(e(k));
    end
    e = sum(sum(e));
end