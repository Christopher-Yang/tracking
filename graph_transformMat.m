function graph_transformMat(data)
groups = {'rot','mir'};
block = {'baseline','pert1','pert2','pert3','pert4','post'};
graph_name = {'Baseline','Early','Train2','Train3','Late','Post'};
Ngroup = length(groups);
Nblock = length(block);
Nsubj = 10;
delay = 50; % number of time steps to delay the target relative to the hand; time step = 1/130.004 sec
% change delay to a vector of values to test for MSE with different delays;
% default delay is 50

paramsInit = [1 0 0 1]; % initialize parameters for optimization

for l = 1:Ngroup
    for i = 1:Nsubj
        for j = 1:Nblock
            dat = data.(groups{l}){i}.(block{j}); % store data into dat
            hand = [dat.Rhand.x_pos dat.Rhand.y_pos]'; % hand position
            target = [dat.target.x_pos dat.target.y_pos]'; % target position
            for k = 1:length(delay) % compute MSE for given delay
                err = @(params) scale(params,hand,target,delay(k)); 
                [params_opt,fval] = fmincon(err,paramsInit); % fit transformation matrix
                transformMat(:,:,k) = [params_opt(1:2); params_opt(3:4)]; % shape into 2x2 matrix
                MSE(k) = fval; % record MSE
            end
            if length(delay) ~= 1 % if testing multiple delays, this chooses delay with lowest MSE
                p = islocalmin(MSE);
                idx(j,i,l) = find(p==1);
            else % if just using one delay, then use that delay
                idx = 1;
            end
            
            % store fitted matrices: the first and second dimensions of
            % each transformMat_* correspond to the transformation matrix
            % for the jth block and the ith subject
            if l == 1
                transformMat_vmr(:,:,j,i) = transformMat(:,:,idx);
            else
                transformMat_mr(:,:,j,i) = transformMat(:,:,idx);
            end
            
            if l == 1
                thetaInit = 0;
                err2 = @(theta) fit_rotMat(theta,transformMat_vmr(:,:,j,i));
                theta_opt(j,i) = fmincon(err2,thetaInit);
            end
        end
    end
end

% estimate rotation angle for the VMR group
for k = 1:Nsubj
    for p = 1:Nblock
        thetaInit = 0;
        err2 = @(theta) fit_rotMat(theta,transformMat_vmr(:,:,j,i));
        theta_opt(j,i) = fmincon(err2,thetaInit);
    end
end

% compute orthogonal gain for mirror-reversal group
R = rotz(-45); % 45 degree clockwise rotation matrix
R = R(1:2,1:2);
perpAxis = R*[1 0]';

for k = 1:Nsubj
    for p = 1:Nblock
        scaleOrth(p,k) = perpAxis'*transformMat_mr(:,:,p,k)*perpAxis;
    end
end

bestDelay = delay(idx).*(1/130.004); % store the delay which minimizes MSE

% This section generates "transformation_matrix.csv" which contains the 
% relevant values of the transformation matrices for statistical analysis 
% in R.
% vmr = cat(3,squeeze(-transformMat_vmr(1,2,[1 5 6],:)),squeeze(transformMat_vmr(2,1,[1 5 6],:)));
% mr = cat(3,squeeze(transformMat_mr(1,2,[1 5 6],:)),squeeze(transformMat_mr(2,1,[1 5 6],:)));
% vmr = reshape(mean(vmr,3)',[30 1]);
% mr = reshape(mean(mr,3)',[30 1]);
% z = [vmr; mr];
% dlmwrite('transformation_matrix.csv',z);

%% generate Figure 3A
% make color maps
col1 = [1 0 0];
col2 = [1 1 1];
Nstep = 100;
map1 = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2),Nstep)', linspace(col1(3),col2(3),Nstep)'];

col1 = [1 1 1];
col2 = [0 0 1];
map2 = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2),Nstep)', linspace(col1(3),col2(3),Nstep)'];

map = [map1; map2];
clims = [-1 1];

% average the transformation matrices across subjects
mat1 = mean(transformMat_vmr,4); % average rotation group
mat2 = mean(transformMat_mr,4); % average mirror-reversal group
gblocks = [1 2 5 6]; % blocks to graph

% plot 2x2 transformation matrices
figure(3); clf
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

% visualize transformation matrices as vectors
for i = 1:4 % compute covariance between target movement in one axis with hand movement in both axes
    x(:,:,i,1) = cov(squeeze(transformMat_vmr(1,1,gblocks(i),:)),squeeze(transformMat_vmr(2,1,gblocks(i),:)));
    y(:,:,i,1) = cov(squeeze(transformMat_vmr(1,2,gblocks(i),:)),squeeze(transformMat_vmr(2,2,gblocks(i),:)));
    x(:,:,i,2) = cov(squeeze(transformMat_mr(1,1,gblocks(i),:)),squeeze(transformMat_mr(2,1,gblocks(i),:)));
    y(:,:,i,2) = cov(squeeze(transformMat_mr(1,2,gblocks(i),:)),squeeze(transformMat_mr(2,2,gblocks(i),:)));
end

figure(4); clf
for i = 1:4
    subplot(2,4,i); hold on
    plot([0 1],[0 0],'k') % unit x vector
    plot([0 0],[0 1],'k') % unit y vector
    a = error_ellipse(x(:,:,i,1),mat1(:,1,gblocks(i)),'color',col1,'conf',0.95); % draw confidence ellipses 
    b = error_ellipse(y(:,:,i,1),mat1(:,2,gblocks(i)),'color',col2,'conf',0.95);
    patch(a.XData,a.YData,col1,'FaceAlpha',0.2) % shade confidence ellipses
    patch(b.XData,b.YData,col2,'FaceAlpha',0.2)
    plot([0 mat1(1,1,gblocks(i))],[0 mat1(2,1,gblocks(i))],'Color',col1,'LineWidth',1.5) % plot vectors
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
    plot([0 1],[0 0],'k') % unit x vector
    plot([0 0],[0 1],'k') % unit y vector
    a = error_ellipse(x(:,:,i,1),mat2(:,1,gblocks(i)),'color',col1,'conf',0.95); % draw confidence ellipses
    b = error_ellipse(y(:,:,i,1),mat2(:,2,gblocks(i)),'color',col2,'conf',0.95);
    patch(a.XData,a.YData,col1,'FaceAlpha',0.2) % shade confidence ellipses
    patch(b.XData,b.YData,col2,'FaceAlpha',0.2)
    plot([0 mat2(1,1,gblocks(i))],[0 mat2(2,1,gblocks(i))],'Color',col1,'LineWidth',1.5) % plot vectors
    plot([0 mat2(1,2,gblocks(i))],[0 mat2(2,2,gblocks(i))],'Color',col2,'LineWidth',1.5)
    axis([-0.9 1.2 -0.9 1.2])
    axis square
    if i == 1
        ylabel('Mirror-Reversal')
    end
    set(gcf,'Renderer','painters')
end

%% generate Figure 3B
col = lines;
col = col(1:7,:);

% average the off-diagonal elements of the transformation matrices
vmr = cat(3,squeeze(-transformMat_vmr(1,2,:,:)),squeeze(transformMat_vmr(2,1,:,:)));
mr = cat(3,squeeze(transformMat_mr(1,2,:,:)),squeeze(transformMat_mr(2,1,:,:)));
vmr = mean(vmr,3);
mr = mean(mr,3);

% variables for plotting
idx = [1 2 0 0 3 4];

figure(5); clf
% plot responses for all tracking blocks
for k = 1:2
    subplot(2,2,k); hold on 
    plot([0 7],[0 0],'--k','LineWidth',1)
    
    % data from individual subjects
    if k == 1
        plot(1:6,vmr,'k','Color',[0 0 0 0.5])
    else
        plot(1:6,mr,'k','Color',[0 0 0 0.5])
    end
    
    % mean across subjects
    for i = 1:6
        if i == 3 || i == 4 % plot blocks between Early and Late in black
            if k == 1
                plot(i,mean(vmr(i,:)),'.k','MarkerSize',24,'LineWidth',1)
            else
                plot(i,mean(mr(i,:)),'.k','MarkerSize',24,'LineWidth',1)
            end
        else % all other blocks plot with colors
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

% plot responses for just baseline and post-learning to assess aftereffects
for k = 3:4
    subplot(2,2,k); hold on
    plot([0 7],[0 0],'--k','LineWidth',1)
    
    % data from individual subjects
    if k == 3
        plot([1 6],vmr([1 6],:),'k','Color',[0 0 0 0.5])
    else
        plot([1 6],mr([1 6],:),'k','Color',[0 0 0 0.5])
    end
    
    % mean across subjects
    for i = [1 6]
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

%% generate Figure 3C

% set variables
gblocks = 1:6; % blocks to plot
thetaMu = mean(theta_opt,2); % mean of rotation group's angle
thetaSE = std(theta_opt,[],2)/sqrt(10); % standard error of angle
scaleOrthMu = mean(scaleOrth,2); % mean of mirror-reversal group's orthogonal scaling
scaleOrthSE = std(scaleOrth,[],2)/sqrt(Nsubj); % standard error of scaling

% set color map
colors = lines;
colors = colors(1:7,:);
idx = [1 2 0 0 3 4];

figure(6); clf
subplot(1,2,1); hold on
plot([0 7],[0 0],'--k','LineWidth',1) % ideal baseline response
plot([0 7],[90 90],'--k','LineWidth',1) % ideal compensation
plot(theta_opt(gblocks,:),'Color',[0 0 0 0.5]) % plot individual subjects
% plot mean across subjects
for i = 1:length(gblocks)
    if i == 3 || i == 4 % plot blocks between Early and Late in black
        plot(i,thetaMu(gblocks(i)),'.','Color',[0 0 0],'MarkerSize',24,'LineWidth',1)
    else % all other blocks plotted in color
        plot(i,thetaMu(gblocks(i)),'.','Color',colors(idx(i),:),'MarkerSize',24,'LineWidth',1)
    end
end
set(gca,'Xcolor','none')
title('Rotation')
xticks([1 2 5 6])
xticklabels(graph_name([1 2 5 6]))
yticks(0:30:90)
ylabel(['Angle (' char(176) ')'])
axis([0.5 6.5 -10 100])

subplot(1,2,2); hold on
plot([0 7],[1 1],'--k','LineWidth',1) % ideal baseline response
plot([0 7],[-1 -1],'--k','LineWidth',1) % ideal compensation
plot([0 7],[0 0],'k','LineWidth',1)
plot(scaleOrth(gblocks,:),'Color',[0 0 0 0.5]) % plot individual subjects
% plot mean across subjects
for i = 1:length(gblocks)
    if i == 3 || i == 4 % plot blocks between Early and Late in black
        plot(i,scaleOrthMu(gblocks(i)),'.','Color',[0 0 0],'MarkerSize',24,'LineWidth',1)
    else% all other blocks plotted in color
        plot(i,scaleOrthMu(gblocks(i)),'.','Color',colors(idx(i),:),'MarkerSize',24,'LineWidth',1)
    end
end
set(gca,'Xcolor','none')
title('Mirror-Reversal')
xticks([1 2 5 6])
xticklabels(graph_name([1 2 5 6]))
yticks(-1:0.5:1)
ylabel('Scaling (orthogonal to mirror axis)')
xlim([0.5 6.5])

%%
function e = scale(params,hand,target,delay) 
    rotMat = [params(1:2); params(3:4)];
    rotTarget = rotMat*target;
    d = (hand(:,delay+1:end)-rotTarget(:,1:end-delay)).^2;
    e = mean(sum(d,1));
end

function e = fit_rotMat(theta,rotMat_opt)
    rot = rotz(theta);
    rot = rot(1:2,1:2);
    e = (rot-rotMat_opt).^2;
    e = sum(sum(e));
end
end