clear all
load dat
groups = {'rot','mir'};
block = {'no_rot1','rot1','rot2','rot3','rot4','no_rot2'};
Nblock = length(block);
Nsubj = 10;
delay = 50;

for l = 1:2
    for i = 1:Nsubj
        for j = 1:Nblock
            dat = data.(groups{l}){i}.(block{j});
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
                rotMat_vmr(:,:,j,i) = rotMat(:,:,idx);
            else
                rotMat_mr(:,:,j,i) = rotMat(:,:,idx);
            end
            
            % estimate rotation angle for the VMR group
            if l == 1
                thetaInit = 0;
                err2 = @(theta) sim_error2(theta,rotMat_vmr(:,:,j,i));
                theta_opt(j,i) = fmincon(err2,thetaInit);
            end
        end
    end
end

bestDelay = delay(idx).*(1/130.004);
thetaBar = mean(theta_opt,2);
thetaSE = std(theta_opt,[],2)/sqrt(10);

% vmr = cat(3,squeeze(-rotMat_vmr(1,2,[1 5 6],:)),squeeze(rotMat_vmr(2,1,[1 5 6],:)));
% mr = cat(3,squeeze(rotMat_mr(1,2,[1 5 6],:)),squeeze(rotMat_mr(2,1,[1 5 6],:)));
% vmr = reshape(mean(vmr,3)',[30 1]);
% mr = reshape(mean(mr,3)',[30 1]);
% z = [vmr; mr];
% dlmwrite('time_matrix.csv',z);
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
mat1 = mean(rotMat_vmr,4);
mat2 = mean(rotMat_mr,4);
gblocks = [1 2 5 6];
figure(1); clf
for i = 1:4
    subplot(2,4,i)
    imagesc(mat1(:,:,gblocks(i)),clims)
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
    imagesc(mat2(:,:,gblocks(i)),clims)
    colormap(map)
    axis square
    if i == 1
        ylabel('Mirror-Reversal')
    end
    set(gca,'Xtick',[],'Ytick',[])
    axis square
end

col1 = [0 128 0]/255;
col2 = [128 0 128]/255;

gblocks = [1 2 5 6];
for i = 1:4
    x(:,:,i,1) = cov(squeeze(rotMat_vmr(1,1,gblocks(i),:)),squeeze(rotMat_vmr(2,1,gblocks(i),:)));
    y(:,:,i,1) = cov(squeeze(rotMat_vmr(1,2,gblocks(i),:)),squeeze(rotMat_vmr(2,2,gblocks(i),:)));
    x(:,:,i,2) = cov(squeeze(rotMat_mr(1,1,gblocks(i),:)),squeeze(rotMat_mr(2,1,gblocks(i),:)));
    y(:,:,i,2) = cov(squeeze(rotMat_mr(1,2,gblocks(i),:)),squeeze(rotMat_mr(2,2,gblocks(i),:)));
end

figure(2); clf
for i = 1:4
    subplot(2,4,i); hold on
    plot([0 1],[0 0],'k')
    plot([0 0],[0 1],'k')
    a = error_ellipse(x(:,:,i,1),mat1(:,1,gblocks(i)),'color',col1,'conf',0.95);
    b = error_ellipse(y(:,:,i,1),mat1(:,2,gblocks(i)),'color',col2,'conf',0.95);
    patch(a.XData,a.YData,col1,'FaceAlpha',0.2)
    patch(b.XData,b.YData,col2,'FaceAlpha',0.2)
    plot([0 mat1(1,1,gblocks(i))],[0 mat1(2,1,gblocks(i))],'Color',col1,'LineWidth',1.5)
    plot([0 mat1(1,2,gblocks(i))],[0 mat1(2,2,gblocks(i))],'Color',col2,'LineWidth',1.5)
    axis([-0.9 1.2 -0.9 1.2])
    axis square
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
    set(gcf,'Renderer','painters')
    
    subplot(2,4,i+4); hold on
    plot([0 1],[0 0],'k')
    plot([0 0],[0 1],'k')
    a = error_ellipse(x(:,:,i,1),mat2(:,1,gblocks(i)),'color',col1,'conf',0.95);
    b = error_ellipse(y(:,:,i,1),mat2(:,2,gblocks(i)),'color',col2,'conf',0.95);
    patch(a.XData,a.YData,col1,'FaceAlpha',0.2)
    patch(b.XData,b.YData,col2,'FaceAlpha',0.2)
    plot([0 mat2(1,1,gblocks(i))],[0 mat2(2,1,gblocks(i))],'Color',col1,'LineWidth',1.5)
    plot([0 mat2(1,2,gblocks(i))],[0 mat2(2,2,gblocks(i))],'Color',col2,'LineWidth',1.5)
    axis([-0.9 1.2 -0.9 1.2])
    axis square
    if i == 1
        ylabel('Mirror-Reversal')
    end
    set(gcf,'Renderer','painters')
end

%% plot aftereffects
col = lines;
col = col(1:7,:);

% vmr = cat(3,squeeze(-rotMat_vmr(1,2,[1 2 5 6],:)),squeeze(rotMat_vmr(2,1,[1 2 5 6],:)));
% mr = cat(3,squeeze(rotMat_mr(1,2,[1 2 5 6],:)),squeeze(rotMat_mr(2,1,[1 2 5 6],:)));
vmr = cat(3,squeeze(-rotMat_vmr(1,2,:,:)),squeeze(rotMat_vmr(2,1,:,:)));
mr = cat(3,squeeze(rotMat_mr(1,2,:,:)),squeeze(rotMat_mr(2,1,:,:)));
vmr = mean(vmr,3);
mr = mean(mr,3);

n = size(vmr,1);
idx = [1 2 0 0 3 4];

figure(3); clf
for k = 1:2
    subplot(2,2,k); hold on
    plot([0 7],[0 0],'--k','LineWidth',1)
    if k == 1
        plot(1:6,vmr,'k','Color',[0 0 0 0.5])
    else
        plot(1:6,mr,'k','Color',[0 0 0 0.5])
    end
    for i = 1:6
        if i == 3 || i == 4
%             if k == 1
%                 errorbar(i,mean(vmr(i,:)),std(vmr(i,:))/sqrt(Nsubj),'.k','MarkerSize',15,'LineWidth',1)
%             else
%                 errorbar(i,mean(mr(i,:)),std(mr(i,:))/sqrt(Nsubj),'.k','MarkerSize',15,'LineWidth',1)
%             end
            if k == 1
                plot(i,mean(vmr(i,:)),'.k','MarkerSize',24,'LineWidth',1)
            else
                plot(i,mean(mr(i,:)),'.k','MarkerSize',24,'LineWidth',1)
            end
        else
%             if k == 1
%                 errorbar(i,mean(vmr(i,:)),std(vmr(i,:))/sqrt(Nsubj),'.','Color',col(idx(i),:),'MarkerSize',15,'LineWidth',1)
%             else
%                 errorbar(i,mean(mr(i,:)),std(mr(i,:))/sqrt(Nsubj),'.','Color',col(idx(i),:),'MarkerSize',15,'LineWidth',1)
%             end
            if k == 1
                plot(i,mean(vmr(i,:)),'.','Color',col(idx(i),:),'MarkerSize',24,'LineWidth',1)
            else
                plot(i,mean(mr(i,:)),'.','Color',col(idx(i),:),'MarkerSize',24,'LineWidth',1)
            end
        end
    end
    axis([0.5 6.5 -.1 0.7])
    set(gca,'Xtick',[],'TickDir','out')
    ylabel('Cross-axis scaling')
    yticks(-1:0.2:1)
end

for k = 3:4
    subplot(2,2,k); hold on
    plot([0 7],[0 0],'--k','LineWidth',1)
    if k == 3
        plot([1 6],vmr([1 6],:),'k','Color',[0 0 0 0.5])
    else
        plot([1 6],mr([1 6],:),'k','Color',[0 0 0 0.5])
    end
    for i = [1 6]
%         if k == 3
%             errorbar(i,mean(vmr(i,:)),std(vmr(i,:))/sqrt(Nsubj),'.','Color',col(idx(i),:),'MarkerSize',15,'LineWidth',1)
%         else
%             errorbar(i,mean(mr(i,:)),std(mr(i,:))/sqrt(Nsubj),'.','Color',col(idx(i),:),'MarkerSize',15,'LineWidth',1)
%         end
        if k == 3
            plot(i,mean(vmr(i,:)),'.','Color',col(idx(i),:),'MarkerSize',24,'LineWidth',1)
        else
            plot(i,mean(mr(i,:)),'.','Color',col(idx(i),:),'MarkerSize',24,'LineWidth',1)
        end
    end
    axis([0.5 6.5 -.1 0.7])
    set(gca,'Xtick',[],'TickDir','out')
    ylabel('Cross-axis scaling')
    yticks(-1:0.2:1)
end

%% plot hand and target trajectories before and after optimization

subj = {'subj14','subj15','subj16','subj19','subj23','subj26','subj27','subj29','subj30','subj34'};
% subj = {'subj17','subj18','subj21','subj22','subj24','subj25','subj28','subj31','subj32','subj33'};
i = 2; % subject number

figure(4); clf
for j = 1:4
    dat = data.rot_i.(subj{i}).(block{j});
    hand = [dat.Rhand.x_pos dat.Rhand.y_pos]';
    target = [dat.target.x_pos dat.target.y_pos]';
    
    rotTarget = rotMat_vmr(:,:,j,i)*target;
    
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

%% plot mirror-reversal compensation perpendicular to mirroring axis
gblocks = 1:6;
colors = lines;
colors = colors(1:7,:);
idx = [1 2 0 0 3 4];

R = rotz(-45);
R = R(1:2,1:2);
perpAxis = R*[1 0]';

for k = 1:Nsubj
    for p = 1:Nblock
        comp(p,k) = perpAxis'*rotMat_mr(:,:,p,k)*perpAxis;
    end
end

compBar = mean(comp,2);
compSE = std(comp,[],2)/sqrt(Nsubj);

figure(3); clf
subplot(2,2,1); hold on
plot([0 7],[0 0],'--k','LineWidth',1)
plot([0 7],[90 90],'--k','LineWidth',1)
plot(theta_opt(gblocks,:),'Color',[0 0 0 0.5])
for i = 1:length(gblocks)
%     if i == 3 || i == 4
%         errorbar(i,thetaBar(gblocks(i)),thetaSE(gblocks(i)),'-ko','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','LineWidth',1)
%     else
%         errorbar(i,thetaBar(gblocks(i)),thetaSE(gblocks(i)),'-o','Color',colors(idx(i),:),'MarkerFaceColor',colors(idx(i),:),'MarkerEdgeColor','none','LineWidth',1)
%     end
    if i == 3 || i == 4
        plot(i,thetaBar(gblocks(i)),'.','Color',[0 0 0],'MarkerSize',24,'LineWidth',1)
    else
        plot(i,thetaBar(gblocks(i)),'.','Color',colors(idx(i),:),'MarkerSize',24,'LineWidth',1)
    end
end
set(gca,'Xcolor','none')
title('Rotation')
xticks(1:6)
xticklabels(graph_name(gblocks))
yticks(0:30:90)
ylabel(['Angle (' char(176) ')'])
axis([0.5 6.5 -10 100])

subplot(2,2,2); hold on
plot([0 7],[1 1],'--k','LineWidth',1)
plot([0 7],[-1 -1],'--k','LineWidth',1)
plot([0 7],[0 0],'k','LineWidth',1)
plot(comp(gblocks,:),'Color',[0 0 0 0.5])
for i = 1:length(gblocks)
%     if i == 3 || i == 4
%         errorbar(i,compBar(gblocks(i)),compSE(gblocks(i)),'-ko','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','LineWidth',1)
%     else
%         errorbar(i,compBar(gblocks(i)),compSE(gblocks(i)),'-o','Color',colors(idx(i),:),'MarkerFaceColor',colors(idx(i),:),'MarkerEdgeColor','none','LineWidth',1)
%     end
    if i == 3 || i == 4
        plot(i,compBar(gblocks(i)),'.','Color',[0 0 0],'MarkerSize',24,'LineWidth',1)
    else
        plot(i,compBar(gblocks(i)),'.','Color',colors(idx(i),:),'MarkerSize',24,'LineWidth',1)
    end
end
set(gca,'Xcolor','none')
title('Mirror-Reversal')
xticks(1:6)
xticklabels(graph_name(gblocks))
yticks(-1:0.5:1)
ylabel('Scaling (orthogonal to mirror axis)')
xlim([0.5 6.5])

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