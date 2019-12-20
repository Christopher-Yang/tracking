clear all
load dat
groups = {'rot','mir'};
blocks = {'no_rot1','rot1','rot2','rot3','rot4','no_rot2'};
names = {'x_x_all','x_y_all','y_x_all','y_y_all'};
Nsubj = length(data.rot)-1;
Nblocks = length(blocks);
for q = 1:2
    for p = 1:Nsubj
        for k = 1:Nblocks
            for i = 1:length(names)
                cplx(:,i,k) = mean(data.(groups{q}){p}.(blocks{k}).phasors.Rhand.(names{i}).ratio,2);
            end
        end
        
        theta = zeros([(Nblocks*4)-1 1]);
        paramsInit = [cplx(:,1,1); theta];
        
        error = @(params) scale(params,cplx); 
        paramsOpt = fmincon(error,paramsInit);
        phasor(:,p,q) = paramsOpt(1:7);
        thetaOpt = [1; paramsOpt(8:end)];
        thetaOpt = reshape(thetaOpt, [4 Nblocks]);
        
        for k = 1:Nblocks
            rotMat(:,:,k,p,q) = reshape(thetaOpt(:,k), [2 2]);
            if q == 1
                thetaInit = 0;
                err2 = @(theta) fit_rotMat(theta,rotMat(:,:,k,p,1));
                thetaFit(k,p) = fmincon(err2,thetaInit);
            end
        end
        
        ph = abs(phasor(:,p,q));
        lambda = mean(ph);
        rotMat(:,:,:,p,q) = rotMat(:,:,:,p,q).*lambda;
        
        if q == 1
            vmrBase(:,:,p) = rotMat(:,:,1,p,q);
            vmrEarly(:,:,p) = rotMat(:,:,2,p,q);
            vmrLate(:,:,p) = rotMat(:,:,5,p,q);
            vmrAfter(:,:,p) = rotMat(:,:,6,p,q);
        else
            mrBase(:,:,p) = rotMat(:,:,1,p,q);
            mrEarly(:,:,p) = rotMat(:,:,2,p,q);
            mrLate(:,:,p) = rotMat(:,:,5,p,q);
            mrAfter(:,:,p) = rotMat(:,:,6,p,q);
        end
    end
end

thetaFit_bar = mean(thetaFit,2);
thetaFit_se = std(thetaFit,[],2)/sqrt(10);

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

%% bar plots of aftereffects
col = lines;
col = col(1:7,:);

vmrXY = [squeeze(vmrBase(1,2,:)) squeeze(vmrLate(1,2,:)) squeeze(vmrAfter(1,2,:))];
vmrYX = [squeeze(vmrBase(2,1,:)) squeeze(vmrLate(2,1,:)) squeeze(vmrAfter(2,1,:))];
mrXY = [squeeze(mrBase(1,2,:)) squeeze(mrLate(1,2,:)) squeeze(mrAfter(1,2,:))];
mrYX = [squeeze(mrBase(2,1,:)) squeeze(mrLate(2,1,:)) squeeze(mrAfter(2,1,:))];


n = size(vmrXY,2);
idx = [1 3 4];

figure(1); clf
subplot(2,2,1); hold on
plot([-1 10],[0 0],'--k','LineWidth',1)
plot(1:3,vmrXY','k','Color',[0 0 0 0.5],'MarkerSize',20)
plot(4:6,mrXY','k','Color',[0 0 0 0.5],'MarkerSize',20)
for i = 1:3
    scatter(repelem(i,Nsubj),vmrXY(:,i),25,col(idx(i),:),'filled','MarkerFaceAlpha',0.5)
    errorbar(i,mean(vmrXY(:,i)),2*std(vmrXY(:,i))/sqrt(Nsubj),'.','Color',col(idx(i),:),'MarkerSize',30,'LineWidth',1)
    
    scatter(repelem(i+n,Nsubj),mrXY(:,i),25,col(idx(i),:),'filled','MarkerFaceAlpha',0.5)
    errorbar(i+n,mean(mrXY(:,i)),2*std(mrXY(:,i))/sqrt(Nsubj),'.','Color',col(idx(i),:),'MarkerSize',30,'LineWidth',1)
end
xlim([0.5 6.5])
ylabel('XY gain')
set(gca,'Xtick',[],'TickDir','out')
yticks(-1:0.2:1)

subplot(2,2,3); hold on
plot([-1 10],[0 0],'--k','LineWidth',1)
plot(1:3,vmrYX','k','Color',[0 0 0 0.5],'MarkerSize',20)
plot(4:6,mrYX','k','Color',[0 0 0 0.5],'MarkerSize',20)
for i = 1:3
    scatter(repelem(i,Nsubj),vmrYX(:,i),25,col(idx(i),:),'filled','MarkerFaceAlpha',0.5)
    errorbar(i,mean(vmrYX(:,i)),2*std(vmrYX(:,i))/sqrt(Nsubj),'.','Color',col(idx(i),:),'MarkerSize',30,'LineWidth',1)
    
    scatter(repelem(i+n,Nsubj),mrYX(:,i),25,col(idx(i),:),'filled','MarkerFaceAlpha',0.5)
    errorbar(i+n,mean(mrYX(:,i)),2*std(mrYX(:,i))/sqrt(Nsubj),'.','Color',col(idx(i),:),'MarkerSize',30,'LineWidth',1)
end
xlim([0.5 6.5])
ylabel('YX gain')
set(gca,'TickDir','out')
xticks([2 5])
yticks(-1:0.2:1)
xticklabels({'Rotation','Mirror-Reversal'})

subplot(2,2,2); hold on
plot([-1 10],[0 0],'--k','LineWidth',1)
plot([1 3],vmrXY(:,[1 3])','k','Color',[0 0 0 0.5],'MarkerSize',20)
plot([4 6],mrXY(:,[1 3])','k','Color',[0 0 0 0.5],'MarkerSize',20)
for i = [1 3]
    scatter(repelem(i,Nsubj),vmrXY(:,i),25,col(idx(i),:),'filled','MarkerFaceAlpha',0.5)
    errorbar(i,mean(vmrXY(:,i)),2*std(vmrXY(:,i))/sqrt(Nsubj),'.','Color',col(idx(i),:),'MarkerSize',30,'LineWidth',1)
    
    scatter(repelem(i+n,Nsubj),mrXY(:,i),25,col(idx(i),:),'filled','MarkerFaceAlpha',0.5)
    errorbar(i+n,mean(mrXY(:,i)),2*std(mrXY(:,i))/sqrt(Nsubj),'.','Color',col(idx(i),:),'MarkerSize',30,'LineWidth',1)
end
xlim([0.5 6.5])
set(gca,'Xtick',[],'TickDir','out')
yticks(-1:0.2:1)

subplot(2,2,4); hold on
plot([-1 10],[0 0],'--k','LineWidth',1)
plot([1 3],vmrYX(:,[1 3])','k','Color',[0 0 0 0.5],'MarkerSize',20)
plot([4 6],mrYX(:,[1 3])','k','Color',[0 0 0 0.5],'MarkerSize',20)
for i = [1 3]
    scatter(repelem(i,Nsubj),vmrYX(:,i),25,col(idx(i),:),'filled','MarkerFaceAlpha',0.5)
    errorbar(i,mean(vmrYX(:,i)),2*std(vmrXY(:,i))/sqrt(Nsubj),'.','Color',col(idx(i),:),'MarkerSize',30,'LineWidth',1)
    
    scatter(repelem(i+n,Nsubj),mrYX(:,i),25,col(idx(i),:),'filled','MarkerFaceAlpha',0.5)
    errorbar(i+n,mean(mrYX(:,i)),2*std(mrYX(:,i))/sqrt(Nsubj),'.','Color',col(idx(i),:),'MarkerSize',30,'LineWidth',1)
end
xlim([0.5 6.5])
set(gca,'TickDir','out')
xticks([2 5])
yticks(-1:0.2:1)
xticklabels({'Rotation','Mirror-Reversal'})

%% plot vectors
% for single subjects
% subj = 5;
% vmrBaseMu = vmrBase(:,:,subj);
% vmrEarlyMu = vmrEarly(:,:,subj);
% vmrLateMu = vmrLate(:,:,subj);
% vmrAfterMu = vmrAfter(:,:,subj);
% 
% mrBaseMu = mrBase(:,:,subj);
% mrEarlyMu = mrEarly(:,:,subj);
% mrLateMu = mrLate(:,:,subj);
% mrAfterMu = mrAfter(:,:,subj);

% for averaging across subjects
% vmrBaseMu = mean(vmrBase,3);
% vmrEarlyMu = mean(vmrEarly,3);
% vmrLateMu = mean(vmrLate,3);
% vmrAfterMu = mean(vmrAfter,3);
% 
% mrBaseMu = mean(mrBase,3);
% mrEarlyMu = mean(mrEarly,3);
% mrLateMu = mean(mrLate,3);
% mrAfterMu = mean(mrAfter,3);

rotMat_mu = squeeze(mean(rotMat,4));

col1 = [0 128 0]/255;
col2 = [128 0 128]/255;

gblocks = [1 2 5 6];
for k = 1:2
    for i = 1:4
        x(:,:,i,k) = cov(squeeze(rotMat(1,1,gblocks(i),:,k)),squeeze(rotMat(2,1,gblocks(i),:,k)));
        y(:,:,i,k) = cov(squeeze(rotMat(1,2,gblocks(i),:,k)),squeeze(rotMat(2,2,gblocks(i),:,k)));
    end
end

figure(1); clf
for k = 1:2
    for i = 1:4
        subplot(2,4,(k-1)*4+i); hold on
        plot([0 1],[0 0],'k')
        plot([0 0],[0 1],'k')
        a = error_ellipse(x(:,:,i,k),rotMat_mu(:,1,gblocks(i),k),'color',col1,'conf',0.95);
        b = error_ellipse(y(:,:,i,k),rotMat_mu(:,2,gblocks(i),k),'color',col2,'conf',0.95);
        patch(a.XData,a.YData,col1,'FaceAlpha',0.2)
        patch(b.XData,b.YData,col2,'FaceAlpha',0.2)
        plot([0 rotMat_mu(1,1,gblocks(i),k)],[0 rotMat_mu(2,1,gblocks(i),k)],'Color',col1,'LineWidth',1.5)
        plot([0 rotMat_mu(1,2,gblocks(i),k)],[0 rotMat_mu(2,2,gblocks(i),k)],'Color',col2,'LineWidth',1.5)
        if k == 1
            if i == 1
                title('Baseline')
                ylabel('Rotation')
            elseif i == 2
                title('Early')
            elseif i == 3
                title('Late')
            elseif i == 4
                title('Post')
            end
        else
            if i == 1
                ylabel('Mirror Reversal')
            end
        end
        axis([-0.75 1.25 -0.75 1.25])
        axis square
        set(gcf,'Renderer','painters')
    end
end

%% fitted phasors for individual subjects
group = 2;
subj = 9;
blockIdx = 6;

plot_subj(data.(groups{group}){subj}.(blocks{blockIdx}).phasors.Rhand)
subplot(2,2,1)
plot(rotMat(1,1,blockIdx,subj,group)*phasor(:,subj,group),'k','LineWidth',3)
xticks(-1:1)
yticks(-1:1)
set(gca,'TickDir','out')
subplot(2,2,2)
plot(rotMat(1,2,blockIdx,subj,group)*phasor(:,subj,group),'k','LineWidth',3)
xticks(-1:1)
yticks(-1:1)
set(gca,'TickDir','out')
subplot(2,2,3)
plot(rotMat(2,1,blockIdx,subj,group)*phasor(:,subj,group),'k','LineWidth',3)
xticks(-1:1)
yticks(-1:1)
set(gca,'TickDir','out')
subplot(2,2,4)
plot(rotMat(2,2,blockIdx,subj,group)*phasor(:,subj,group),'k','LineWidth',3)
xticks(-1:1)
yticks(-1:1)
set(gca,'TickDir','out')

%% gain matrices
figure(1); clf
subplot(2,4,1)
imagesc(mean(vmrBase,3),clims)
colormap(map)
set(gca,'Xtick',[],'Ytick',[])
axis square
title('Baseline')
ylabel('Rotation')

subplot(2,4,2)
imagesc(mean(vmrEarly,3),clims)
colormap(map)
set(gca,'Xtick',[],'Ytick',[])
axis square
title('Early')

subplot(2,4,3)
imagesc(mean(vmrLate,3),clims)
colormap(map)
set(gca,'Xtick',[],'Ytick',[])
axis square
title('Late')

subplot(2,4,4)
imagesc(mean(vmrAfter,3),clims)
colormap(map)
set(gca,'Xtick',[],'Ytick',[])
axis square
title('After')

subplot(2,4,5)
imagesc(mean(mrBase,3),clims)
colormap(map)
set(gca,'TickDir','out','Xtick',[],'Ytick',[])
axis square
ylabel('Mirror Reversal')

subplot(2,4,6)
imagesc(mean(mrEarly,3),clims)
% imagesc([1 0; 0 1],clims)
colormap(map)
set(gca,'TickDir','out','Xtick',[],'Ytick',[])
axis square

subplot(2,4,7)
imagesc(mean(mrLate,3),clims)
% imagesc([0 1; -1 0],clims)
colormap(map)
set(gca,'Xtick',[],'Ytick',[])
axis square

subplot(2,4,8)
imagesc(mean(mrAfter,3),clims)
% imagesc([0 1; 1 0],clims)
colormap(map)
set(gca,'Xtick',[],'Ytick',[])
axis square

%%
function e = scale(params,cplx_data)
    cplx_fit = params(1:7);
    theta = [1; params(8:end)];
    theta = reshape(theta,[4 length(theta)/4]);
    for k = 1:size(cplx_data,3)
        rotMat = reshape(theta(:,k),[2 2]);
        for i = 1:4
            gain = rotMat(i);
            phasor = gain*cplx_fit;
            e(i,k) = norm(phasor - cplx_data(:,i,k))^2;
        end
    end
    e = sum(sum(e));
end

function e = fit_rotMat(theta,rotMat_opt)
    rot = rotz(theta);
    rot = rot(1:2,1:2);
    e = (rot-rotMat_opt).^2;
    e = sum(sum(e));
end