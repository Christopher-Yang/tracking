Nsubj = length(data)-1;
blocks = {'baseline','dark_baseline','rot1','rot2','rot3','rot4','after'};
names = {'x_x_all','x_y_all','y_x_all','y_y_all'};
output = 'Rhand';

for p = 1:Nsubj
    for k = 1:length(blocks)
        for i = 1:length(names)
            cplx(:,i,k) = mean(data{p}.(blocks{k}).phasors.(output).(names{i}).ratio,2);
        end
    end
    
    theta = zeros([27 1]);
    paramsInit = [cplx(:,1,1); theta];
    
    error = @(params) scale(params,cplx);
    paramsOpt = fmincon(error,paramsInit);
    phasor(:,p) = paramsOpt(1:7);
    thetaOpt = [1; paramsOpt(8:end)];
    thetaOpt = reshape(thetaOpt, [4 length(blocks)]);
    
    for k = 1:length(blocks)
        rotMat(:,:,k,p) = reshape(thetaOpt(:,k), [2 2]);
        thetaInit = 0;
        err2 = @(theta) fit_rotMat(theta,rotMat(:,:,k,p,1));
        thetaFit(k,p) = fmincon(err2,thetaInit);
    end
    
    ph = abs(phasor(:,p));
    lambda = mean(ph);
    
    Base(:,:,p) = rotMat(:,:,1,p).*lambda;
    Dark_base(:,:,p) = rotMat(:,:,2,p).*lambda;
    Early(:,:,p) = rotMat(:,:,3,p).*lambda;
    Late(:,:,p) = rotMat(:,:,6,p).*lambda;
    After(:,:,p) = rotMat(:,:,7,p).*lambda;
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

vmrXY = [squeeze(Base(1,2,:)) squeeze(Late(1,2,:)) squeeze(After(1,2,:))];
vmrYX = [squeeze(Base(2,1,:)) squeeze(Late(2,1,:)) squeeze(After(2,1,:))];

figure(1); clf
subplot(2,4,2:3); 
hold on
plot([-1 10],[0 0],'--k','LineWidth',1)
plot(1:0.5:2,vmrXY','k','Color',[0 0 0 0.5],'MarkerSize',20)
scatter(repelem(1,Nsubj),vmrXY(:,1),25,col(1,:),'filled','MarkerFaceAlpha',0.5)
scatter(repelem(1.5,Nsubj),vmrXY(:,2),25,col(3,:),'filled','MarkerFaceAlpha',0.5)
scatter(repelem(2,Nsubj),vmrXY(:,3),25,col(4,:),'filled','MarkerFaceAlpha',0.5)
errorbar(1,mean(vmrXY(:,1)),2*std(vmrXY(:,1))/sqrt(10),'.','Color',col(1,:),'MarkerSize',30,'LineWidth',1)
errorbar(1.5,mean(vmrXY(:,2)),2*std(vmrXY(:,2))/sqrt(10),'.','Color',col(3,:),'MarkerSize',30,'LineWidth',1)
errorbar(2,mean(vmrXY(:,3)),2*std(vmrXY(:,3))/sqrt(10),'.','Color',col(4,:),'MarkerSize',30,'LineWidth',1)

xlim([0.5 2.5])
ylabel('XY gain')
set(gca,'Xtick',[],'TickDir','out')
yticks(-1:0.2:1)

subplot(2,4,6:7); 
hold on
plot([-1 10],[0 0],'--k','LineWidth',1)
plot(1:0.5:2,vmrYX','k','Color',[0 0 0 0.5],'MarkerSize',20)
scatter(repelem(1,Nsubj),vmrYX(:,1),25,col(1,:),'filled','MarkerFaceAlpha',0.5)
scatter(repelem(1.5,Nsubj),vmrYX(:,2),25,col(3,:),'filled','MarkerFaceAlpha',0.5)
scatter(repelem(2,Nsubj),vmrYX(:,3),25,col(4,:),'filled','MarkerFaceAlpha',0.5)
errorbar(1,mean(vmrYX(:,1)),2*std(vmrYX(:,1))/sqrt(10),'.','Color',col(1,:),'MarkerSize',30,'LineWidth',1)
errorbar(1.5,mean(vmrYX(:,2)),2*std(vmrYX(:,2))/sqrt(10),'.','Color',col(3,:),'MarkerSize',30,'LineWidth',1)
errorbar(2,mean(vmrYX(:,3)),2*std(vmrYX(:,3))/sqrt(10),'.','Color',col(4,:),'MarkerSize',30,'LineWidth',1)

xlim([0.5 2.5])
ylabel('YX gain')
set(gca,'TickDir','out')
yticks(-1:0.2:1)

%% plot vectors
% for single subjects
subj = 1;
BaseMu = Base(:,:,subj);
DarkBaseMu = Dark_base(:,:,subj);
EarlyMu = Early(:,:,subj);
LateMu = Late(:,:,subj);
AfterMu = After(:,:,subj);

% for averaging across subjects
% BaseMu = mean(Base,3);
% DarkBaseMu = mean(Dark_base,3);
% EarlyMu = mean(Early,3);
% LateMu = mean(Late,3);
% AfterMu = mean(After,3);

% gain matrices
figure(1); clf
subplot(2,5,1)
imagesc(BaseMu,clims)
colormap(map)
set(gca,'Xtick',[],'Ytick',[])
axis square
title('Baseline')

subplot(2,5,2)
imagesc(DarkBaseMu,clims)
colormap(map)
set(gca,'Xtick',[],'Ytick',[])
axis square
title('Dark Baseline')

subplot(2,5,3)
imagesc(EarlyMu,clims)
colormap(map)
set(gca,'Xtick',[],'Ytick',[])
axis square
title('Early')

subplot(2,5,4)
imagesc(LateMu,clims)
colormap(map)
set(gca,'Xtick',[],'Ytick',[])
axis square
title('Late')

subplot(2,5,5)
imagesc(AfterMu,clims)
colormap(map)
set(gca,'Xtick',[],'Ytick',[])
axis square
title('After')

% vectors
subplot(2,5,6); hold on
plot([0 BaseMu(1,1)],[0 BaseMu(2,1)],'LineWidth',1.5)
plot([0 BaseMu(1,2)],[0 BaseMu(2,2)],'LineWidth',1.5)
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
title('Baseline')
axis([-0.45 1 -0.45 1])
axis square

subplot(2,5,7); hold on
plot([0 DarkBaseMu(1,1)],[0 DarkBaseMu(2,1)],'LineWidth',1.5)
plot([0 DarkBaseMu(1,2)],[0 DarkBaseMu(2,2)],'LineWidth',1.5)
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
title('Dark Baseline')
axis([-0.45 1 -0.45 1])
axis square

subplot(2,5,8); hold on
plot([0 EarlyMu(1,1)],[0 EarlyMu(2,1)],'LineWidth',1.5)
plot([0 EarlyMu(1,2)],[0 EarlyMu(2,2)],'LineWidth',1.5)
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
title('Early')
axis([-0.45 1 -0.45 1])
axis square

subplot(2,5,9); hold on
plot([0 LateMu(1,1)],[0 LateMu(2,1)],'LineWidth',1.5)
plot([0 LateMu(1,2)],[0 LateMu(2,2)],'LineWidth',1.5)
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
title('Late')
axis([-0.45 1 -0.45 1])
axis square

subplot(2,5,10); hold on
plot([0 AfterMu(1,1)],[0 AfterMu(2,1)],'LineWidth',1.5)
plot([0 AfterMu(1,2)],[0 AfterMu(2,2)],'LineWidth',1.5)
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
title('Post')
axis([-0.45 1 -0.45 1])
axis square

%% fitted phasors for individual subjects
subj = 2;
blockIdx = 1;

plot_subj(data{subj}.(blocks{blockIdx}).phasors.(output),2)
subplot(2,2,1)
plot(rotMat(1,1,blockIdx,subj)*phasor(:,subj),'k','LineWidth',3)
subplot(2,2,2)
plot(rotMat(1,2,blockIdx,subj)*phasor(:,subj),'k','LineWidth',3)
subplot(2,2,3)
plot(rotMat(2,1,blockIdx,subj)*phasor(:,subj),'k','LineWidth',3)
subplot(2,2,4)
plot(rotMat(2,2,blockIdx,subj)*phasor(:,subj),'k','LineWidth',3)

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