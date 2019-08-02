clear all
load dat
subj_rot = {'subj17','subj18','subj21','subj22','subj24','subj25','subj28','subj31','subj32','subj33'};
subj_rot_i = {'subj14','subj15','subj16','subj19','subj23','subj26','subj27','subj29','subj30','subj34'};
blocks = {'no_rot1','rot1','rot4','no_rot2'};
names = {'x_x_all','x_y_all','y_x_all','y_y_all'};
for q = 1:2
    for p = 1:length(subj_rot)
        for k = 1:length(blocks)
            for i = 1:length(names)
                if q == 1
                    cplx(:,i,k) = mean(data.rot.(subj_rot{p}).(blocks{k}).(names{i}).ratio,2);
                else
                    cplx(:,i,k) = mean(data.rot_i.(subj_rot_i{p}).(blocks{k}).(names{i}).ratio,2);
                end
            end
        end
        
        theta = zeros([15 1]);
        paramsInit = [cplx(:,1,1); theta];
        
        error = @(params) scale(params,cplx); 
        paramsOpt = fmincon(error,paramsInit);
        phasor(:,p,q) = paramsOpt(1:7);
        thetaOpt = [1; paramsOpt(8:end)];
        thetaOpt = reshape(thetaOpt, [4 4]);
        
        for k = 1:length(thetaOpt)
            rotMat(:,:,k,p,q) = reshape(thetaOpt(:,k), [2 2]);
        end
        
        ph = abs(phasor(:,p,q));
        lambda = mean(ph);
        
        if q == 1
            vmrBase(:,:,p) = rotMat(:,:,1,p,q).*lambda;
            vmrEarly(:,:,p) = rotMat(:,:,2,p,q).*lambda;
            vmrLate(:,:,p) = rotMat(:,:,3,p,q).*lambda;
            vmrAfter(:,:,p) = rotMat(:,:,4,p,q).*lambda;
        else
            mrBase(:,:,p) = rotMat(:,:,1,p,q).*lambda;
            mrEarly(:,:,p) = rotMat(:,:,2,p,q).*lambda;
            mrLate(:,:,p) = rotMat(:,:,3,p,q).*lambda;
            mrAfter(:,:,p) = rotMat(:,:,4,p,q).*lambda;
        end
    end
end

%%
col = lines;
col = col(1:7,:);

vmrXY = [squeeze(vmrBase(1,2,:)) squeeze(vmrAfter(1,2,:))];
vmrYX = [squeeze(vmrBase(2,1,:)) squeeze(vmrAfter(2,1,:))];
mrXY = [squeeze(mrBase(1,2,:)) squeeze(mrAfter(1,2,:))];
mrYX = [squeeze(mrBase(2,1,:)) squeeze(mrAfter(2,1,:))];

figure(1); clf
subplot(2,3,2); 
hold on
plot([-1 10],[0 0],'--k','LineWidth',1)
plot(1:2,vmrXY','k','Color',[0 0 0 0.5],'MarkerSize',20)
scatter(repelem(1,10),vmrXY(:,1),25,col(1,:),'filled','MarkerFaceAlpha',0.5)
scatter(repelem(2,10),vmrXY(:,2),25,col(4,:),'filled','MarkerFaceAlpha',0.5)
errorbar(1,mean(vmrXY(:,1)),2*std(vmrXY(:,1))/sqrt(10),'.','Color',col(1,:),'MarkerSize',30,'LineWidth',1)
errorbar(2,mean(vmrXY(:,2)),2*std(vmrXY(:,2))/sqrt(10),'.','Color',col(4,:),'MarkerSize',30,'LineWidth',1)

plot(3.5:4.5,mrXY','k','Color',[0 0 0 0.5],'MarkerSize',20)
scatter(repelem(3.5,10),mrXY(:,1),25,col(1,:),'filled','MarkerFaceAlpha',0.5)
scatter(repelem(4.5,10),mrXY(:,2),25,col(4,:),'filled','MarkerFaceAlpha',0.5)
errorbar(3.5,mean(mrXY(:,1)),2*std(mrXY(:,1))/sqrt(10),'.','Color',col(1,:),'MarkerSize',30,'LineWidth',1)
errorbar(4.5,mean(mrXY(:,2)),2*std(mrXY(:,2))/sqrt(10),'.','Color',col(4,:),'MarkerSize',30,'LineWidth',1)

xlim([0.5 5])
ylabel('XY gain')
set(gca,'Xtick',[],'TickDir','out')

subplot(2,3,5); 
hold on
plot([-1 10],[0 0],'--k','LineWidth',1)
plot(1:2,vmrYX','k','Color',[0 0 0 0.5],'MarkerSize',20)
scatter(repelem(1,10),vmrYX(:,1),25,col(1,:),'filled','MarkerFaceAlpha',0.5)
scatter(repelem(2,10),vmrYX(:,2),25,col(4,:),'filled','MarkerFaceAlpha',0.5)
errorbar(1,mean(vmrYX(:,1)),2*std(vmrYX(:,1))/sqrt(10),'.','Color',col(1,:),'MarkerSize',30,'LineWidth',1)
errorbar(2,mean(vmrYX(:,2)),2*std(vmrYX(:,2))/sqrt(10),'.','Color',col(4,:),'MarkerSize',30,'LineWidth',1)

plot(3.5:4.5,mrYX','k','Color',[0 0 0 0.5],'MarkerSize',20)
scatter(repelem(3.5,10),mrYX(:,1),25,col(1,:),'filled','MarkerFaceAlpha',0.5)
scatter(repelem(4.5,10),mrYX(:,2),25,col(4,:),'filled','MarkerFaceAlpha',0.5)
errorbar(3.5,mean(mrYX(:,1)),2*std(mrYX(:,1))/sqrt(10),'.','Color',col(1,:),'MarkerSize',30,'LineWidth',1)
errorbar(4.5,mean(mrYX(:,2)),2*std(mrYX(:,2))/sqrt(10),'.','Color',col(4,:),'MarkerSize',30,'LineWidth',1)

xlim([0.5 5])
ylabel('YX gain')
set(gca,'TickDir','out')
xticks([1.5 4])
xticklabels({'VMR','MR'})

%%
vmrBaseMu = mean(vmrBase,3);
vmrEarlyMu = mean(vmrEarly,3);
vmrLateMu = mean(vmrLate,3);
vmrAfterMu = mean(vmrAfter,3);

mrBaseMu = mean(mrBase,3);
mrEarlyMu = mean(mrEarly,3);
mrLateMu = mean(mrLate,3);
mrAfterMu = mean(mrAfter,3);

figure(1); clf
subplot(2,4,1); hold on
plot([0 vmrBaseMu(1,1)],[0 vmrBaseMu(2,1)],'LineWidth',1.5)
plot([0 vmrBaseMu(1,2)],[0 vmrBaseMu(2,2)],'LineWidth',1.5)
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
title('Baseline')
ylabel('Rotation')
axis([-0.45 1 -0.45 1])
axis square

subplot(2,4,2); hold on
plot([0 vmrEarlyMu(1,1)],[0 vmrEarlyMu(2,1)],'LineWidth',1.5)
plot([0 vmrEarlyMu(1,2)],[0 vmrEarlyMu(2,2)],'LineWidth',1.5)
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
title('Early')
axis([-0.45 1 -0.45 1])
axis square

subplot(2,4,3); hold on
plot([0 vmrLateMu(1,1)],[0 vmrLateMu(2,1)],'LineWidth',1.5)
plot([0 vmrLateMu(1,2)],[0 vmrLateMu(2,2)],'LineWidth',1.5)
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
title('Late')
axis([-0.45 1 -0.45 1])
axis square

subplot(2,4,4); hold on
plot([0 vmrAfterMu(1,1)],[0 vmrAfterMu(2,1)],'LineWidth',1.5)
plot([0 vmrAfterMu(1,2)],[0 vmrAfterMu(2,2)],'LineWidth',1.5)
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
title('Post')
axis([-0.45 1 -0.45 1])
axis square

subplot(2,4,5); hold on
plot([0 mrBaseMu(1,1)],[0 mrBaseMu(2,1)],'LineWidth',1.5)
plot([0 mrBaseMu(1,2)],[0 mrBaseMu(2,2)],'LineWidth',1.5)
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
ylabel('Mirror Reversal')
axis([-0.45 1 -0.45 1])
axis square

subplot(2,4,6); hold on
plot([0 mrEarlyMu(1,1)],[0 mrEarlyMu(2,1)],'LineWidth',1.5)
plot([0 mrEarlyMu(1,2)],[0 mrEarlyMu(2,2)],'LineWidth',1.5)
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
axis([-0.45 1 -0.45 1])
axis square

subplot(2,4,7); hold on
plot([0 mrLateMu(1,1)],[0 mrLateMu(2,1)],'LineWidth',1.5)
plot([0 mrLateMu(1,2)],[0 mrLateMu(2,2)],'LineWidth',1.5)
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
axis([-0.45 1 -0.45 1])
axis square

subplot(2,4,8); hold on
plot([0 mrAfterMu(1,1)],[0 mrAfterMu(2,1)],'LineWidth',1.5)
plot([0 mrAfterMu(1,2)],[0 mrAfterMu(2,2)],'LineWidth',1.5)
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
axis([-0.45 1 -0.45 1])
axis square

%%
group = 1;
subj = 4;
blockIdx = 4;

plot_subj(data.rot.(subj_rot{subj}).(blocks{blockIdx}),0)
subplot(2,2,1)
plot(rotMat(1,1,blockIdx,subj,group)*phasor(:,subj,group),'k','LineWidth',3)
subplot(2,2,2)
plot(rotMat(1,2,blockIdx,subj,group)*phasor(:,subj,group),'k','LineWidth',3)
subplot(2,2,3)
plot(rotMat(2,1,blockIdx,subj,group)*phasor(:,subj,group),'k','LineWidth',3)
subplot(2,2,4)
plot(rotMat(2,2,blockIdx,subj,group)*phasor(:,subj,group),'k','LineWidth',3)

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
imagesc(mean(vmrBase,3),clims)
colormap(map)
set(gca,'TickDir','out','Xtick',[],'Ytick',[])
axis square
title('Baseline')
ylabel('Rotation')

subplot(2,4,2)
imagesc(mean(vmrEarly,3),clims)
colormap(map)
set(gca,'TickDir','out','Xtick',[],'Ytick',[])
axis square
title('Early')

subplot(2,4,3)
imagesc(mean(vmrLate,3),clims)
colormap(map)
set(gca,'TickDir','out','Xtick',[],'Ytick',[])
axis square
title('Late')

subplot(2,4,4)
imagesc(mean(vmrAfter,3),clims)
colormap(map)
set(gca,'TickDir','out','Xtick',[],'Ytick',[])
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