% clear all
% load dat
groups = {'rot','rot_i'};
block = {'no_rot1','rot1','rot4','no_rot2'};
delay = 50;

for l = 1:2
    if l == 1
        subj = {'subj17','subj18','subj21','subj22','subj24','subj25','subj28','subj31','subj32','subj33'};
    else
        subj = {'subj14','subj15','subj16','subj19','subj23','subj26','subj27','subj29','subj30','subj34'};
    end
    for i = 1:10
        for j = 1:4
            dat = data.(groups{l}).(subj{i}).(block{j});
            hand = [dat.Rhand.x_pos dat.Rhand.y_pos]';
            target = [dat.target.x_pos dat.target.y_pos]';
            paramsInit = [1 0 0 1];
            for k = 1:length(delay)
                err = @(params) sim_error(params,hand,target,delay(k));
                [params_opt,fval] = fmincon(err,paramsInit);
                rotMat(:,:,k) = [params_opt(1:2); params_opt(3:4)];
                MSE(k) = fval;
            end
            if length(delay) ~= 1
                p = islocalmin(MSE);
                idx(j,i,l) = find(p==1);
            else
                idx = 1;
            end
            
            % store fitted matrices
            if l == 1
                rotMat_opt_vmr(:,:,j,i) = rotMat(:,:,idx);
            else
                rotMat_opt_mr(:,:,j,i) = rotMat(:,:,idx);
            end
            
            % estimate rotation angle for the VMR group
            if l == 1
                thetaInit = 0;
                err2 = @(theta) sim_error2(theta,rotMat_opt_vmr(:,:,j,i));
                theta_opt(j,i) = fmincon(err2,thetaInit);
            end
        end
    end
end

bestDelay = delay(idx).*(1/130.004);
theta_bar = mean(theta_opt,2);
theta_se = std(theta_opt,[],2)/sqrt(10);

% vmr12 = reshape(squeeze(rotMat_opt_vmr(1,2,[1 4],:))',[20 1]);
% vmr21 = reshape(squeeze(rotMat_opt_vmr(2,1,[1 4],:))',[20 1]);
% mr12 = reshape(squeeze(rotMat_opt_mr(1,2,[1 4],:))',[20 1]);
% mr21 = reshape(squeeze(rotMat_opt_mr(2,1,[1 4],:))',[20 1]);
% z = [vmr12; vmr21; mr12; mr21];
% dlmwrite('time_matrix2.csv',z);
%% plot fitted matrices as well as column vector representation
col1 = [1 0 0];
col2 = [1 1 1];
Nstep = 100;
map1 = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2),Nstep)', linspace(col1(3),col2(3),Nstep)'];

col1 = [1 1 1];
col2 = [0 0 1];
map2 = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2),Nstep)', linspace(col1(3),col2(3),Nstep)'];

map = [map1; map2];
clims = [-1 1];

% subj = 10;
% mat1 = rotMat_opt_vmr(:,:,:,subj);
% mat2 = rotMat_opt_mr(:,:,:,subj);
mat1 = mean(rotMat_opt_vmr,4);
mat2 = mean(rotMat_opt_mr,4);
figure(1); clf
for i = 1:4
    subplot(2,4,i)
    imagesc(mat1(:,:,i),clims)
    colormap(map)
    axis square
    if i == 1
        title('Baseline')
        ylabel('Rotation')
    elseif i == 2
        title('Early')
    elseif i == 3
        title('Late')
    else
        title('After')
    end
    set(gca,'Xtick',[],'Ytick',[])
    axis square
    
    subplot(2,4,i+4)
    imagesc(mat2(:,:,i),clims)
    colormap(map)
    axis square
    if i == 1
        ylabel('Mirror-Reversal')
    end
    set(gca,'Xtick',[],'Ytick',[])
    axis square
end

figure(2); clf
for i = 1:4
    subplot(2,4,i); hold on
    plot([0 mat1(1,1,i)],[0 mat1(2,1,i)],'LineWidth',1.5)
    plot([0 mat1(1,2,i)],[0 mat1(2,2,i)],'LineWidth',1.5)
    plot([0 1],[0 0],'k')
    plot([0 0],[0 1],'k')
    axis([-0.65 1 -0.65 1])
    axis square
    if i == 1
        ylabel('Rotation')
    end
    
    subplot(2,4,i+4); hold on
    plot([0 mat2(1,1,i)],[0 mat2(2,1,i)],'LineWidth',1.5)
    plot([0 mat2(1,2,i)],[0 mat2(2,2,i)],'LineWidth',1.5)
    plot([0 1],[0 0],'k')
    plot([0 0],[0 1],'k')
    axis([-0.65 1 -0.65 1])
    axis square
    if i == 1
        ylabel('Mirror-Reversal')
    end
end

%% plot aftereffects
col = lines;
col = col(1:7,:);

vmr12 = squeeze(rotMat_opt_vmr(1,2,[1 3 4],:));
vmr21 = squeeze(rotMat_opt_vmr(2,1,[1 3 4],:));
mr12 = squeeze(rotMat_opt_mr(1,2,[1 3 4],:));
mr21 = squeeze(rotMat_opt_mr(2,1,[1 3 4],:));

figure(3); clf
subplot(2,4,1:2); hold on
plot([-1 10],[0 0],'--k','LineWidth',1)
plot(1:0.5:2,vmr12,'k','Color',[0 0 0 0.5],'MarkerSize',20)
scatter(repelem(1,10),vmr12(1,:),25,col(1,:),'filled','MarkerFaceAlpha',0.5)
scatter(repelem(1.5,10),vmr12(2,:),25,col(3,:),'filled','MarkerFaceAlpha',0.5)
scatter(repelem(2,10),vmr12(3,:),25,col(4,:),'filled','MarkerFaceAlpha',0.5)
errorbar(1,mean(vmr12(1,:)),2*std(vmr12(1,:))/sqrt(10),'.','Color',col(1,:),'MarkerSize',30,'LineWidth',1)
errorbar(1.5,mean(vmr12(2,:)),2*std(vmr12(2,:))/sqrt(10),'.','Color',col(3,:),'MarkerSize',30,'LineWidth',1)
errorbar(2,mean(vmr12(3,:)),2*std(vmr12(3,:))/sqrt(10),'.','Color',col(4,:),'MarkerSize',30,'LineWidth',1)

plot(3:0.5:4,mr12,'k','Color',[0 0 0 0.5],'MarkerSize',20)
scatter(repelem(3,10),mr12(1,:),25,col(1,:),'filled','MarkerFaceAlpha',0.5)
scatter(repelem(3.5,10),mr12(2,:),25,col(3,:),'filled','MarkerFaceAlpha',0.5)
scatter(repelem(4,10),mr12(3,:),25,col(4,:),'filled','MarkerFaceAlpha',0.5)
errorbar(3,mean(mr12(1,:)),2*std(mr12(1,:))/sqrt(10),'.','Color',col(1,:),'MarkerSize',30,'LineWidth',1)
errorbar(3.5,mean(mr12(2,:)),2*std(mr12(2,:))/sqrt(10),'.','Color',col(3,:),'MarkerSize',30,'LineWidth',1)
errorbar(4,mean(mr12(3,:)),2*std(mr12(3,:))/sqrt(10),'.','Color',col(4,:),'MarkerSize',30,'LineWidth',1)
xlim([0.5 4.5])
ylabel('Element 1,2')
set(gca,'Xtick',[],'TickDir','out')
yticks(-1:0.2:1)

subplot(2,4,3:4); hold on
plot([-1 10],[0 0],'--k','LineWidth',1)
plot(1:2,vmr12([1 3],:),'k','Color',[0 0 0 0.5],'MarkerSize',20)
scatter(repelem(1,10),vmr12(1,:),25,col(1,:),'filled','MarkerFaceAlpha',0.5)
scatter(repelem(2,10),vmr12(3,:),25,col(4,:),'filled','MarkerFaceAlpha',0.5)
errorbar(1,mean(vmr12(1,:)),2*std(vmr12(1,:))/sqrt(10),'.','Color',col(1,:),'MarkerSize',30,'LineWidth',1)
errorbar(2,mean(vmr12(3,:)),2*std(vmr12(3,:))/sqrt(10),'.','Color',col(4,:),'MarkerSize',30,'LineWidth',1)

plot(3:4,mr12([1 3],:),'k','Color',[0 0 0 0.5],'MarkerSize',20)
scatter(repelem(3,10),mr12(1,:),25,col(1,:),'filled','MarkerFaceAlpha',0.5)
scatter(repelem(4,10),mr12(3,:),25,col(4,:),'filled','MarkerFaceAlpha',0.5)
errorbar(3,mean(mr12(1,:)),2*std(mr12(1,:))/sqrt(10),'.','Color',col(1,:),'MarkerSize',30,'LineWidth',1)
errorbar(4,mean(mr12(3,:)),2*std(mr12(3,:))/sqrt(10),'.','Color',col(4,:),'MarkerSize',30,'LineWidth',1)
xlim([0.5 4.5])
set(gca,'Xtick',[],'TickDir','out')
yticks(-1:0.2:1)

subplot(2,4,5:6); hold on
plot([-1 10],[0 0],'--k','LineWidth',1)
plot(1:0.5:2,vmr21,'k','Color',[0 0 0 0.5],'MarkerSize',20)
scatter(repelem(1,10),vmr21(1,:),25,col(1,:),'filled','MarkerFaceAlpha',0.5)
scatter(repelem(1.5,10),vmr21(2,:),25,col(3,:),'filled','MarkerFaceAlpha',0.5)
scatter(repelem(2,10),vmr21(3,:),25,col(4,:),'filled','MarkerFaceAlpha',0.5)
errorbar(1,mean(vmr21(1,:)),2*std(vmr21(1,:))/sqrt(10),'.','Color',col(1,:),'MarkerSize',30,'LineWidth',1)
errorbar(1.5,mean(vmr21(2,:)),2*std(vmr21(2,:))/sqrt(10),'.','Color',col(3,:),'MarkerSize',30,'LineWidth',1)
errorbar(2,mean(vmr21(3,:)),2*std(vmr21(3,:))/sqrt(10),'.','Color',col(4,:),'MarkerSize',30,'LineWidth',1)

plot(3:0.5:4,mr21,'k','Color',[0 0 0 0.5],'MarkerSize',20)
scatter(repelem(3,10),mr21(1,:),25,col(1,:),'filled','MarkerFaceAlpha',0.5)
scatter(repelem(3.5,10),mr21(2,:),25,col(3,:),'filled','MarkerFaceAlpha',0.5)
scatter(repelem(4,10),mr21(3,:),25,col(4,:),'filled','MarkerFaceAlpha',0.5)
errorbar(3,mean(mr21(1,:)),2*std(mr21(1,:))/sqrt(10),'.','Color',col(1,:),'MarkerSize',30,'LineWidth',1)
errorbar(3.5,mean(mr21(2,:)),2*std(mr21(2,:))/sqrt(10),'.','Color',col(3,:),'MarkerSize',30,'LineWidth',1)
errorbar(4,mean(mr21(3,:)),2*std(mr21(3,:))/sqrt(10),'.','Color',col(4,:),'MarkerSize',30,'LineWidth',1)
axis([0.5 4.5 -0.2 0.65])
ylabel('Element 2,1')
set(gca,'Xtick',[],'TickDir','out')
xticks([1.5 3.5])
yticks(-1:0.2:1)
xticklabels({'Rotation','Mirror-Reversal'})

subplot(2,4,7:8); hold on
plot([-1 10],[0 0],'--k','LineWidth',1)
plot(1:2,vmr21([1 3],:),'k','Color',[0 0 0 0.5],'MarkerSize',20)
scatter(repelem(1,10),vmr21(1,:),25,col(1,:),'filled','MarkerFaceAlpha',0.5)
scatter(repelem(2,10),vmr21(3,:),25,col(4,:),'filled','MarkerFaceAlpha',0.5)
errorbar(1,mean(vmr21(1,:)),2*std(vmr21(1,:))/sqrt(10),'.','Color',col(1,:),'MarkerSize',30,'LineWidth',1)
errorbar(2,mean(vmr21(3,:)),2*std(vmr21(3,:))/sqrt(10),'.','Color',col(4,:),'MarkerSize',30,'LineWidth',1)

plot(3:4,mr21([1 3],:),'k','Color',[0 0 0 0.5],'MarkerSize',20)
scatter(repelem(3,10),mr21(1,:),25,col(1,:),'filled','MarkerFaceAlpha',0.5)
scatter(repelem(4,10),mr21(3,:),25,col(4,:),'filled','MarkerFaceAlpha',0.5)
errorbar(3,mean(mr21(1,:)),2*std(mr21(1,:))/sqrt(10),'.','Color',col(1,:),'MarkerSize',30,'LineWidth',1)
errorbar(4,mean(mr21(3,:)),2*std(mr21(3,:))/sqrt(10),'.','Color',col(4,:),'MarkerSize',30,'LineWidth',1)
xlim([0.5 4.5])
set(gca,'Xtick',[],'TickDir','out')
xticks([1.5 3.5])
yticks(-1:0.2:1)
xticklabels({'Rotation','Mirror-Reversal'})

%% plot hand and target trajectories before and after optimization

subj = {'subj14','subj15','subj16','subj19','subj23','subj26','subj27','subj29','subj30','subj34'};
% subj = {'subj17','subj18','subj21','subj22','subj24','subj25','subj28','subj31','subj32','subj33'};
i = 2; % subject number

figure(4); clf
for j = 1:4
    dat = data.rot_i.(subj{i}).(block{j});
    hand = [dat.Rhand.x_pos dat.Rhand.y_pos]';
    target = [dat.target.x_pos dat.target.y_pos]';
    
    rotTarget = rotMat_opt_vmr(:,:,j,i)*target;
    
    subplot(2,4,j); hold on
    plot(target(1,:),target(2,:))
    plot(hand(1,:),hand(2,:))
    axis square
    if j == 1
        ylabel('Before Optimization')
    end
    
    subplot(2,4,j+4); hold on
    plot(rotTarget(1,:),rotTarget(2,:))
    plot(hand(1,:),hand(2,:))
    axis square
    if j == 1
        ylabel('After Optimization')
    end
end
%%
function e = sim_error(params,hand,target,delay)
    rotMat = [params(1:2); params(3:4)];
    rotTarget = rotMat*target;
    d = (hand(:,delay+1:end)-rotTarget(:,1:end-delay)).^2;
    e = mean(sum(d,1));
end

function e = sim_error2(theta,rotMat_opt)
    rot = rotz(theta);
    rot = rot(1:2,1:2);
    e = (rot-rotMat_opt).^2;
    e = sum(sum(e));
end