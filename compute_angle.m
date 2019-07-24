% clear all
% load dat
subj_rot = {'subj11','subj12','subj13','subj14','subj15','subj17','subj18','subj20','subj21','subj22'};
blocks = {'baseline','rot1','rot4','after'};
names = {'x_x_all','x_y_all','y_x_all','y_y_all'};

for p = 1:length(subj_rot)
    for k = 1:length(blocks)
        for i = 1:length(names)
            cplx(:,i,k) = mean(data.rot.(subj_rot{p}).(blocks{k}).(names{i}).ratio,2);
        end
    end
    
    theta = zeros([15 1]);
    paramsInit = [cplx(:,1,1); theta];
    
    error = @(params) scale(params,cplx);
    paramsOpt = fmincon(error,paramsInit);
    phasor(:,p) = paramsOpt(1:7);
    thetaOpt = [1; paramsOpt(8:end)];
    thetaOpt = reshape(thetaOpt, [4 4]);
    
    for k = 1:length(thetaOpt)
        rotMat(:,:,k,p) = reshape(thetaOpt(:,k), [2 2]);
    end
    
    ph = abs(phasor(:,p));
    lambda = mean(ph);
    
    base(:,:,p) = rotMat(:,:,1,p).*lambda;
    early(:,:,p) = rotMat(:,:,2,p).*lambda;
    late(:,:,p) = rotMat(:,:,3,p).*lambda;
    after(:,:,p) = rotMat(:,:,4,p).*lambda;
end

%% 
col = lines;
col = col(1:7,:);

xy = [squeeze(base(1,2,:)) squeeze(after(1,2,:))];
yx = [squeeze(base(2,1,:)) squeeze(after(2,1,:))];

figure(1); clf
hold on
plot([-1 10],[0 0],'--k','LineWidth',1)
plot(1:2,xy','k','Color',[0 0 0 0.5],'MarkerSize',20)
scatter(repelem(1,10),xy(:,1),25,col(1,:),'filled','MarkerFaceAlpha',0.5)
scatter(repelem(2,10),xy(:,2),25,col(4,:),'filled','MarkerFaceAlpha',0.5)
errorbar(1,mean(xy(:,1)),2*std(xy(:,1))/sqrt(10),'.','Color',col(1,:),'MarkerSize',30,'LineWidth',1)
errorbar(2,mean(xy(:,2)),2*std(xy(:,2))/sqrt(10),'.','Color',col(4,:),'MarkerSize',30,'LineWidth',1)

plot(3.5:4.5,yx','k','Color',[0 0 0 0.5],'MarkerSize',20)
scatter(repelem(3.5,10),yx(:,1),25,col(1,:),'filled','MarkerFaceAlpha',0.5)
scatter(repelem(4.5,10),yx(:,2),25,col(4,:),'filled','MarkerFaceAlpha',0.5)
errorbar(3.5,mean(yx(:,1)),2*std(yx(:,1))/sqrt(10),'.','Color',col(1,:),'MarkerSize',30,'LineWidth',1)
errorbar(4.5,mean(yx(:,2)),2*std(yx(:,2))/sqrt(10),'.','Color',col(4,:),'MarkerSize',30,'LineWidth',1)
ylabel('Gain')
set(gca,'Xtick',[],'TickDir','out')
xlim([0 5.5])
xticks([1.5 4])
xticklabels({'XY','YX'})

%%
baseMu = mean(base,3);
earlyMu = mean(early,3);
lateMu = mean(late,3);
afterMu = mean(after,3);

figure(1); clf
subplot(1,4,1); hold on
plot([0 baseMu(1,1)],[0 baseMu(2,1)],'LineWidth',1.5)
plot([0 baseMu(1,2)],[0 baseMu(2,2)],'LineWidth',1.5)
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
title('Baseline')
ylabel('Rotation')
axis([-0.45 1 -0.45 1])
axis square

subplot(1,4,2); hold on
plot([0 earlyMu(1,1)],[0 earlyMu(2,1)],'LineWidth',1.5)
plot([0 earlyMu(1,2)],[0 earlyMu(2,2)],'LineWidth',1.5)
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
title('Early')
axis([-0.45 1 -0.45 1])
axis square

subplot(1,4,3); hold on
plot([0 lateMu(1,1)],[0 lateMu(2,1)],'LineWidth',1.5)
plot([0 lateMu(1,2)],[0 lateMu(2,2)],'LineWidth',1.5)
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
title('Late')
axis([-0.45 1 -0.45 1])
axis square

subplot(1,4,4); hold on
plot([0 afterMu(1,1)],[0 afterMu(2,1)],'LineWidth',1.5)
plot([0 afterMu(1,2)],[0 afterMu(2,2)],'LineWidth',1.5)
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
title('Post')
axis([-0.45 1 -0.45 1])
axis square

%%
subj = 10;
blockIdx = 4;

plot_subj(data.rot.(subj_rot{subj}).(blocks{blockIdx}),0)
subplot(2,2,1)
plot(rotMat(1,1,blockIdx,subj)*phasor(:,subj),'k','LineWidth',3)
subplot(2,2,2)
plot(rotMat(1,2,blockIdx,subj)*phasor(:,subj),'k','LineWidth',3)
subplot(2,2,3)
plot(rotMat(2,1,blockIdx,subj)*phasor(:,subj),'k','LineWidth',3)
subplot(2,2,4)
plot(rotMat(2,2,blockIdx,subj)*phasor(:,subj),'k','LineWidth',3)

%%
col1 = [1 0 0];
col2 = [1 1 1];
Nstep = 100;
map1 = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2),Nstep)', linspace(col1(3),col2(3),Nstep)'];

col1 = [1 1 1];
col2 = [0 0 1];
map2 = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2),Nstep)', linspace(col1(3),col2(3),Nstep)'];

map = [map1; map2];
clims = [-1 1];

figure(1); clf
subplot(2,4,1)
imagesc(mean(base,3),clims)
colormap(map)
set(gca,'TickDir','out','Xtick',[],'Ytick',[])
axis square
title('Baseline')
ylabel('Rotation')

subplot(2,4,2)
imagesc(mean(early,3),clims)
colormap(map)
set(gca,'TickDir','out','Xtick',[],'Ytick',[])
axis square
title('Early')

subplot(2,4,3)
imagesc(mean(late,3),clims)
colormap(map)
set(gca,'TickDir','out','Xtick',[],'Ytick',[])
axis square
title('Late')

subplot(2,4,4)
imagesc(mean(after,3),clims)
colormap(map)
set(gca,'TickDir','out','Xtick',[],'Ytick',[])
axis square
title('After')

%%
function e = scale(params,cplx_data)
    cplx_fit = params(1:7);
    theta = [1; params(8:end)];
    theta = reshape(theta,[4 4]);
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