function graph_transformMat(data)
% performs analysis and makes plots for transformation matrices

% set variables for analysis
groups = {'rot','mir'};
block = {'baseline','pert1','pert2','pert3','pert4','post'};
graph_name = {'Baseline','Early','Train2','Train3','Late','Post'};
Ngroup = length(groups);
Nblock = length(block);
Nsubj = 10;
Ntrials = 8;

% delay.mat contains the time delay that minimizes the mean-squared error
% between target and hand trajectories for each trial. Loading this data 
% file considerably reduces computation time. This was calculated using
% lines ......... in graph_transformMat.m. To generate delay.mat, comment
% lines 18-19, uncomment lines 24-25, and run graph_transformMat.m.
load delay
delayTest = 0;

% Number of time steps to delay the target relative to the hand (each time 
% step is 1/130.004 sec; default delay is 50). Change delay to a vector of 
% values to test the mean-squared error with different delays.
% delay = 10:10:200; 
% delayTest = 1

paramsInit = [1 0 0 1]; % initialize parameters for optimization
transformMat_vmr = NaN(2,2,Ntrials,Nblock,Nsubj); % preallocate matrices
transformMat_mr = NaN(2,2,Ntrials,Nblock,Nsubj);

for l = 1:Ngroup
    for i = 1:Nsubj
        for j = 1:Nblock
            dat = data.(groups{l}){i}.(block{j});
            for m = 1:Ntrials
                % hand position
                hand = [dat.Rhand.x_pos(:,m) dat.Rhand.y_pos(:,m)]';
                % target position
                target = [dat.target.x_pos(:,m) dat.target.y_pos(:,m)]';
                
                if delayTest == 1
                    transformMat = NaN(2,2,Ntrials,Nblock);
                    for k = 1:length(delay) % compute MSE for given delay
                        err = @(params) scale(params,hand,target,delay(k));
                        
                        % fit transformation matrix
                        [params_opt,fval] = fmincon(err,paramsInit);
                        
                        % shape into 2x2 matrix
                        transformMat(:,:,m,k) = [params_opt(1:2); params_opt(3:4)];
                        
                        % record MSE
                        MSE(k,m,j,i,l) = fval;
                    end
                    
                    % select the transformation matrix that minimizes MSE
                    p = islocalmin(MSE);
                    idx(j,i,l) = find(p==1);
                    transformMat = transformMat(:,:,m,idx);
                else
                    err = @(params) scale(params,hand,target,delay(m,j,i,l));
                    
                    % fit transformation matrix
                    params_opt = fmincon(err,paramsInit);
                    
                    % shape into 2x2 matrix
                    transformMat = [params_opt(1:2); params_opt(3:4)];
                end
                
                % store fitted matrices: the first and second dimensions of
                % each transformMat_* correspond to the transformation matrix
                % for the jth block and the ith subject
                if l == 1
                    transformMat_vmr(:,:,m,j,i) = transformMat;
                else
                    transformMat_mr(:,:,m,j,i) = transformMat;
                end
            end
        end
    end
end

% estimate rotation angle for the VMR group
for i = 1:Nsubj
    for j = 1:Nblock
        for m = 1:Ntrials
            H = transformMat_vmr(:,:,m,j,i)';
            [U,S,V] = svd(H);
            if det(H') >= 0
                R = V*U';
                theta_opt(m,j,i) = atan2(R(2,1),R(1,1))*180/pi;
            else
                theta_opt(m,j,i) = NaN;
            end
        end
    end
end

% compute orthogonal gain for mirror-reversal group
R = rotz(-45); % 45 degree clockwise rotation matrix
R = R(1:2,1:2);
perpAxis = R*[1 0]';

for i = 1:Nsubj
    for j = 1:Nblock
        for m = 1:Ntrials
            scaleOrth(m,j,i) = perpAxis'*transformMat_mr(:,:,m,j,i)*perpAxis;
        end
    end
end

%% generate Figure 3A
% make color maps
col1 = [1 0 0];
col2 = [1 1 1];
Nstep = 100;
map1 = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2),...
    Nstep)', linspace(col1(3),col2(3),Nstep)'];

col1 = [1 1 1];
col2 = [0 0 1];
map2 = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2),...
    Nstep)', linspace(col1(3),col2(3),Nstep)'];

map = [map1; map2];
clims = [-1 1];

% average the transformation matrices across subjects
mat1 = mean(transformMat_vmr,5); % average rotation group
mat2 = mean(transformMat_mr,5); % average mirror-reversal group
gblocks = [1 2 5 6]; % blocks to graph

% plot 2x2 transformation matrices
figure(3); clf
for i = 1:4 % iterate over blocks
    for j = 1:Ntrials % iterate over trials
        subplot(4,Ntrials,(i-1)*Ntrials+j)
        imagesc(mat1(:,:,j,gblocks(i)),clims)
        colormap(map)
        axis square
        if j == 1
            if i == 1
                ylabel('Baseline')
            elseif i == 2
                ylabel('Early')
            elseif i == 3
                ylabel('Late')
            else
                ylabel('After')
            end
        end
        set(gca,'Xtick',[],'Ytick',[])
        axis square
    end
end

figure(4); clf
for i = 1:4 % iterate over blocks
    for j = 1:Ntrials % iterate over trials
        subplot(4,Ntrials,(i-1)*Ntrials+j)
        imagesc(mat2(:,:,j,gblocks(i)),clims)
        colormap(map)
        axis square
        if j == 1
            if i == 1
                ylabel('Baseline')
            elseif i == 2
                ylabel('Early')
            elseif i == 3
                ylabel('Late')
            else
                ylabel('After')
            end
        end
        set(gca,'Xtick',[],'Ytick',[])
        axis square
    end
end

%%

col1 = [0 128 0]/255;
col2 = [128 0 128]/255;

% compute covariance between target movement in one axis with hand movement
% in both axes
for i = 1:4
    for j = 1:Ntrials
        x(:,:,j,i,1) = cov(squeeze(transformMat_vmr(1,1,j,gblocks(i),:)),...
            squeeze(transformMat_vmr(2,1,j,gblocks(i),:)));
        y(:,:,j,i,1) = cov(squeeze(transformMat_vmr(1,2,j,gblocks(i),:)),...
            squeeze(transformMat_vmr(2,2,j,gblocks(i),:)));
        x(:,:,j,i,2) = cov(squeeze(transformMat_mr(1,1,j,gblocks(i),:)),...
            squeeze(transformMat_mr(2,1,j,gblocks(i),:)));
        y(:,:,j,i,2) = cov(squeeze(transformMat_mr(1,2,j,gblocks(i),:)),...
            squeeze(transformMat_mr(2,2,j,gblocks(i),:)));
    end
end

% visualize transformation matrices as vectors
figure(5); clf
for i = 1:4
    for j = 1:Ntrials
        subplot(4,Ntrials,(i-1)*Ntrials+j); hold on
        plot([0 1],[0 0],'k') % unit x vector
        plot([0 0],[0 1],'k') % unit y vector
        a = error_ellipse(x(:,:,j,i,1),mat1(:,1,j,gblocks(i)),'color',col1,...
            'conf',0.95); % draw confidence ellipses
        b = error_ellipse(y(:,:,j,i,1),mat1(:,2,j,gblocks(i)),'color',col2,...
            'conf',0.95);
        patch(a.XData,a.YData,col1,'FaceAlpha',0.2) % shade confidence ellipses
        patch(b.XData,b.YData,col2,'FaceAlpha',0.2)
        plot([0 mat1(1,1,j,gblocks(i))],[0 mat1(2,1,j,gblocks(i))],'Color',col1,...
            'LineWidth',1.5) % plot vectors
        plot([0 mat1(1,2,j,gblocks(i))],[0 mat1(2,2,j,gblocks(i))],'Color',col2,...
            'LineWidth',1.5)
        axis([-1 1.3 -1 1.3])
        axis square
        if j == 1
            if i == 1
                ylabel('Baseline')
            elseif i == 2
                ylabel('Early')
            elseif i == 3
                ylabel('Late')
            elseif i == 4
                ylabel('Post')
            end
        end
        set(gcf,'Renderer','painters')
    end
end

figure(6); clf
for i = 1:4
    for j = 1:Ntrials
        subplot(4,Ntrials,(i-1)*Ntrials+j); hold on
        plot([0 1],[0 0],'k') % unit x vector
        plot([0 0],[0 1],'k') % unit y vector
        a = error_ellipse(x(:,:,j,i,2),mat2(:,1,j,gblocks(i)),'color',col1,...
            'conf',0.95); % draw confidence ellipses
        b = error_ellipse(y(:,:,j,i,2),mat2(:,2,j,gblocks(i)),'color',col2,...
            'conf',0.95);
        patch(a.XData,a.YData,col1,'FaceAlpha',0.2) % shade confidence ellipses
        patch(b.XData,b.YData,col2,'FaceAlpha',0.2)
        plot([0 mat2(1,1,j,gblocks(i))],[0 mat2(2,1,j,gblocks(i))],'Color',col1,...
            'LineWidth',1.5) % plot vectors
        plot([0 mat2(1,2,j,gblocks(i))],[0 mat2(2,2,j,gblocks(i))],'Color',col2,...
            'LineWidth',1.5)
        axis([-1 1.3 -1 1.3])
        axis square
        if j == 1
            if i == 1
                ylabel('Baseline')
            elseif i == 2
                ylabel('Early')
            elseif i == 3
                ylabel('Late')
            elseif i == 4
                ylabel('Post')
            end
        end
        set(gcf,'Renderer','painters')
    end
end

%% generate Figure 3B
col = lines;
col = col(1:7,:);

% average the off-diagonal elements of the transformation matrices
vmr = cat(4,squeeze(-transformMat_vmr(1,2,:,:,:)),squeeze...
    (transformMat_vmr(2,1,:,:,:)));
mr = cat(4,squeeze(transformMat_mr(1,2,:,:,:)),squeeze...
    (transformMat_mr(2,1,:,:,:)));
vmr = mean(vmr,4);
mr = mean(mr,4);

% variables for plotting
colIdx = [1 2 0 0 3 4];
vmr = reshape(vmr,[Ntrials*Nblock Nsubj]);
mr = reshape(mr,[Ntrials*Nblock Nsubj]);

% This section generates "transformation_matrix.csv" which contains the 
% relevant values of the transformation matrices for statistical analysis 
% in R.
% z = [vmr(1,:)'; vmr(40,:)'; vmr(41,:)'; mr(1,:)'; mr(40,:)'; mr(41,:)'];
% dlmwrite('transformation_matrix.csv',z);

figure(7); clf
% plot responses for all tracking blocks
for k = 1:2
    subplot(2,2,k); hold on 
    rectangle('Position',[9 -10 31 110],'FaceColor',[0 0 0 0.1],'EdgeColor','none');
    plot([0 Ntrials*Nblock+1],[0 0],'k')
    
    % mean across subjects
    for i = 1:6
        idx = Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials;
        if i == 3 || i == 4 % plot blocks between Early and Late in black
            if k == 1
                plot(idx,vmr(idx,:),'k','Color',[0 0 0 0.5])
                plot(idx,mean(vmr(idx,:),2),'k','LineWidth',3)
            else
                plot(idx,mr(idx,:),'k','Color',[0 0 0 0.5])
                plot(idx,mean(mr(idx,:),2),'k','LineWidth',3)
            end
        else % all other blocks plot with colors
            if k == 1
                plot(idx,vmr(idx,:),'k','Color',[0 0 0 0.5])
                plot(idx,mean(vmr(idx,:),2),'Color',col(colIdx(i),:)...
                    ,'LineWidth',3)
            else
                plot(idx,mr(idx,:),'k','Color',[0 0 0 0.5])
                plot(idx,mean(mr(idx,:),2),'Color',col(colIdx(i),:)...
                    ,'LineWidth',3)
            end
        end
    end
    axis([0.5 Ntrials*Nblock+.5 -.2 0.8])
    set(gca,'TickDir','out')
    xticks(1:8:41)
    ylabel('Off-diagonal scaling')
    yticks(0:0.4:0.8)
end

% plot responses for just baseline and post-learning to assess aftereffects
for k = 3:4
    subplot(2,2,k); hold on
    plot([0 Ntrials*Nblock+1],[0 0],'k')
    idx = [1:8 41:48];
    
    % data from individual subjects
    if k == 3
        plot(idx,vmr(idx,:),'k','Color',[0 0 0 0.5])
    else
        plot(idx,mr(idx,:),'k','Color',[0 0 0 0.5])
    end
    
    % mean across subjects
    for i = [1 6]
        idx = Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials;
        
        if k == 3
            plot(idx,mean(vmr(idx,:),2),'Color',col(colIdx(i),:)...
                ,'LineWidth',3)
        else
            plot(idx,mean(mr(idx,:),2),'Color',col(colIdx(i),:)...
                ,'LineWidth',3)
        end
    end
    axis([0.5 Ntrials*Nblock+0.5 -.2 0.8])
    set(gca,'TickDir','out')
    xticks(1:8:41)
    xlabel('Trial Number')
    ylabel('Off-diagonal scaling')
    yticks(0:0.4:0.8)
end

%% generate Figure 3C

theta_opt2 = reshape(theta_opt,[Nblock*Ntrials Nsubj]);
scaleOrth2 = reshape(scaleOrth,[Nblock*Ntrials Nsubj]);

% set variables
thetaMu = nanmean(theta_opt2,2); % mean of rotation group's angle
scaleOrthMu = mean(scaleOrth2,2); % mean of MR group's orthogonal scaling

% set color map
colors = lines;
colors = colors(1:7,:);
idx = [1 2 0 0 3 4];

figure(8); clf
subplot(1,2,1); hold on
rectangle('Position',[9 -10 31 110],'FaceColor',[0 0 0 0.1],'EdgeColor','none');
plot([0 Nblock*Ntrials+1],[0 0],'--k','LineWidth',1) % ideal baseline response
plot([0 Nblock*Ntrials+1],[90 90],'--k','LineWidth',1) % ideal compensation
% plot mean across subjects
for i = 1:Nblock
    plotIdx = Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials;
    plot(plotIdx,theta_opt2(plotIdx,:),'Color',[0 0 0 0.5],'LineWidth',0.2) % plot individual subjects
    if i == 3 || i == 4 % plot blocks between Early and Late in black
        plot(plotIdx,thetaMu(plotIdx),'Color',[0 0 0],'LineWidth',1.5)
    else % all other blocks plotted in color
        plot(plotIdx,thetaMu(plotIdx),'Color',colors(idx(i),:),'LineWidth',1.5)
    end
end
set(gca,'TickDir','out')
title('Rotation')
xticks(1:8:41)
xlabel('Trial Number')
yticks(0:30:90)
ylabel(['Angle (' char(176) ')'])
axis([0.5 length(thetaMu) + 0.5 -10 100])

subplot(1,2,2); hold on
rectangle('Position',[9 -10 31 110],'FaceColor',[0 0 0 0.1],'EdgeColor','none');
plot([0 Nblock*Ntrials+1],[1 1],'--k','LineWidth',1) % ideal baseline response
plot([0 Nblock*Ntrials+1],[-1 -1],'--k','LineWidth',1) % ideal compensation
plot([0 Nblock*Ntrials+1],[0 0],'k')
% plot mean across subjects
for i = 1:Nblock
    plotIdx = Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials;
    plot(plotIdx,scaleOrth2(plotIdx,:),'Color',[0 0 0 0.5],'LineWidth',0.2) % plot individual subjects
    if i == 3 || i == 4 % plot blocks between Early and Late in black
        plot(plotIdx,scaleOrthMu(plotIdx),'Color',[0 0 0],'LineWidth',1.5)
    else% all other blocks plotted in color
        plot(plotIdx,scaleOrthMu(plotIdx),'Color',colors(idx(i),:),'LineWidth',1.5)
    end
end
set(gca,'TickDir','out')
title('Mirror-Reversal')
xticks(1:8:41)
xlabel('Trial Number')
% xticks([1 2 5 6])
% xticklabels(graph_name([1 2 5 6]))
yticks(-1:0.5:1)
ylabel('Scaling (orthogonal to mirror axis)')
axis([0.5 length(scaleOrthMu)+0.5 -1.1 1.1])

%%
function e = scale(params,hand,target,delay) 
    rotMat = [params(1:2); params(3:4)];
    rotTarget = rotMat*target;
    d = (hand(:,delay+1:end)-rotTarget(:,1:end-delay)).^2;
    e = mean(sum(d,1));
end
end