% Performs analysis and makes plots for alignment matrices. The optional 
% argument "delay_opt" contains the time delay that minimizes the mean-
% squared error between target and hand trajectories for each trial. If
% delay_opt is provided, it considerably reduces computation time. If
% delay_opt is not provided, it will be computed from scratch. You can save
% a new delay_opt by uncommenting line 109.
% 
% This function also generates "alignment_matrix.csv," which contains the
% off-diagonal values of the alignment matrices for statistical analysis in
% R. To generate this file, uncomment lines 315-316.

function graph_alignMatrix(data, delay_opt)

% set variables for analysis
groups = {'rot','mir'}; % names of groups
block = {'baseline','pert1','pert2','pert3','pert4','post'}; % block names
Ngroup = length(groups); % number of groups of subjects
Nblock = length(block); % number of blocks
Nsubj = 10; % number of subjects
Ntrials = 8; % number of trials per block
paramsInit = [1 0 0 1]; % initialize parameters for optimization
alignMat_vmr = NaN(2,2,Ntrials,Nblock,Nsubj); % 
alignMat_mr = NaN(2,2,Ntrials,Nblock,Nsubj);

% This if statement determines whether to use the provided delay_opt or
% compute it from scratch. If computing from scratch, the code tests
% the time delays contained in "delay." This variable is in time steps of
% 1/130 sec. The code then selects which delay minimizes the mean squared
% error between the target and hand.
if nargin > 1
    delayTest = 0;
else
    delay = 10:10:280; 
    delayTest = 1;
end

% this loop computes the alignment matrices
for n = 1:Ngroup % loop over groups
    for i = 1:Nsubj % loop over subjects
        for j = 1:Nblock % loop over blocks
            dat = data.(groups{n}){i}.(block{j});
            for m = 1:Ntrials % loop over trials
                % hand position
                hand = [dat.Rhand.x_pos(:,m) dat.Rhand.y_pos(:,m)]';
                
                % target position
                target = [dat.target.x_pos(:,m) dat.target.y_pos(:,m)]';
                
                % tests various delays if delay_opt is not provided
                if delayTest == 1
                    alignMat = NaN(2,2,Ntrials,Nblock);
                    for k = 1:length(delay) % compute MSE for given delay
                        err = @(params) scale(params,hand,target,delay(k));
                        
                        % fit alignment matrix
                        [params_opt,fval] = fmincon(err,paramsInit);
                        
                        % shape into 2x2 matrix
                        alignMat(:,:,m,k) = [params_opt(1:2); ...
                            params_opt(3:4)];
                        
                        % record MSE
                        MSE(k) = fval;
                    end
                    
                    % select the alignment matrix that minimizes MSE
                    p = islocalmin(MSE);
                    mins = find(p==1);
                    
                    if length(mins) == 0
                        delay_opt(m,j,i,n) = NaN;
                        alignMat = NaN(2);
                    elseif length(mins) == 1
                        idx = mins;
                        delay_opt(m,j,i,n) = delay(mins);
                        alignMat = alignMat(:,:,m,mins);
                    else
                        [~, id] = min(MSE(p));
                        idx = mins(id);
                        delay_opt(m,j,i,n) = delay(idx);
                        alignMat = alignMat(:,:,m,idx);
                    end
                    
                % if delay_opt is provided, uses this for analysis
                else
                    err = @(params) scale(params,hand,target,...
                        delay_opt(m,j,i,n));
                    
                    % fit alignment matrix
                    params_opt = fmincon(err,paramsInit);
                    
                    % shape into 2x2 matrix
                    alignMat = [params_opt(1:2); params_opt(3:4)];
                end
                
                % store fitted matrices: the first and second dimensions of
                % each alignMat_* correspond to the alignment matrix for 
                % the jth block and the ith subject
                if n == 1
                    alignMat_vmr(:,:,m,j,i) = alignMat;
                else
                    alignMat_mr(:,:,m,j,i) = alignMat;
                end
            end
        end
    end
end

% Save delay_opt as a .m file. Uncomment the line below if you would like
% to save a new delay_opt file
% save delay_opt delay_opt

% estimate rotation angle for the VMR group
for i = 1:Nsubj
    for j = 1:Nblock
        for m = 1:Ntrials
            H = alignMat_vmr(:,:,m,j,i)';
            [U,S,V] = svd(H);
            if det(H') >= 0 % if determinant >= 0, rotation matrix can be computed
                R = V*U';
                theta_opt(m,j,i) = atan2(R(2,1),R(1,1))*180/pi;
            else % if rotation matrix can't be computed, set value to NaN
                theta_opt(m,j,i) = NaN;
            end
        end
    end
end

% compute orthogonal gain for mirror-reversal group
R = rotz(-45); % 45 degree clockwise rotation matrix
R = R(1:2,1:2);
perpAxis = R*[1 0]'; % axis perpendicular to the mirroring axis

for i = 1:Nsubj
    for j = 1:Nblock
        for m = 1:Ntrials
            scaleOrth(m,j,i) = perpAxis'*alignMat_mr(:,:,m,j,i)*perpAxis;
        end
    end
end

%% generate Figure 3A
% make color maps for plots
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

% average the alignment matrices across subjects
mat1 = mean(alignMat_vmr,5); % average rotation group
mat2 = mean(alignMat_mr,5); % average mirror-reversal group
gblocks = [1 2 5 6]; % blocks to graph

% plot 2x2 alignment matrices for rotation group
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

% plot 2x2 alignment matrices for mirror-reversal group
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

col1 = [0 128 0]/255;
col2 = [128 0 128]/255;

% compute covariance between target movement in one axis with hand movement
% in both axes
for i = 1:4
    for j = 1:Ntrials
        x(:,:,j,i,1) = cov(squeeze(alignMat_vmr(1,1,j,gblocks(i),:)),...
            squeeze(alignMat_vmr(2,1,j,gblocks(i),:)));
        y(:,:,j,i,1) = cov(squeeze(alignMat_vmr(1,2,j,gblocks(i),:)),...
            squeeze(alignMat_vmr(2,2,j,gblocks(i),:)));
        x(:,:,j,i,2) = cov(squeeze(alignMat_mr(1,1,j,gblocks(i),:)),...
            squeeze(alignMat_mr(2,1,j,gblocks(i),:)));
        y(:,:,j,i,2) = cov(squeeze(alignMat_mr(1,2,j,gblocks(i),:)),...
            squeeze(alignMat_mr(2,2,j,gblocks(i),:)));
    end
end

% visualize alignment matrices as vectors for rotation group
figure(5); clf
for i = 1:4
    for j = 1:Ntrials
        subplot(4,Ntrials,(i-1)*Ntrials+j); hold on
        plot([0 1],[0 0],'k') % unit x vector
        plot([0 0],[0 1],'k') % unit y vector
        a = error_ellipse(x(:,:,j,i,1),mat1(:,1,j,gblocks(i)),'color',...
            col1,'conf',0.95); % draw confidence ellipses
        b = error_ellipse(y(:,:,j,i,1),mat1(:,2,j,gblocks(i)),'color',...
            col2,'conf',0.95);
        patch(a.XData,a.YData,col1,'FaceAlpha',0.2) % shade confidence ellipses
        patch(b.XData,b.YData,col2,'FaceAlpha',0.2)
        plot([0 mat1(1,1,j,gblocks(i))],[0 mat1(2,1,j,gblocks(i))],...
            'Color',col1,'LineWidth',1.5) % plot vectors
        plot([0 mat1(1,2,j,gblocks(i))],[0 mat1(2,2,j,gblocks(i))],...
            'Color',col2,'LineWidth',1.5)
        axis([-1.5 1.4 -1.5 1.4])
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

% visualize alignment matrices as vectors for mirror-reversal group
figure(6); clf
for i = 1:4
    for j = 1:Ntrials
        subplot(4,Ntrials,(i-1)*Ntrials+j); hold on
        plot([0 1],[0 0],'k') % unit x vector
        plot([0 0],[0 1],'k') % unit y vector
        a = error_ellipse(x(:,:,j,i,2),mat2(:,1,j,gblocks(i)),'color',...
            col1,'conf',0.95); % draw confidence ellipses
        b = error_ellipse(y(:,:,j,i,2),mat2(:,2,j,gblocks(i)),'color',...
            col2,'conf',0.95);
        patch(a.XData,a.YData,col1,'FaceAlpha',0.2) % shade confidence ellipses
        patch(b.XData,b.YData,col2,'FaceAlpha',0.2)
        plot([0 mat2(1,1,j,gblocks(i))],[0 mat2(2,1,j,gblocks(i))],...
            'Color',col1,'LineWidth',1.5) % plot vectors
        plot([0 mat2(1,2,j,gblocks(i))],[0 mat2(2,2,j,gblocks(i))],...
            'Color',col2,'LineWidth',1.5)
        axis([-1.5 1.4 -1.5 1.4])
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

% average the off-diagonal elements of the alignment matrices
vmr = cat(4,squeeze(-alignMat_vmr(1,2,:,:,:)),squeeze...
    (alignMat_vmr(2,1,:,:,:)));
mr = cat(4,squeeze(alignMat_mr(1,2,:,:,:)),squeeze...
    (alignMat_mr(2,1,:,:,:)));
vmr = mean(vmr,4);
mr = mean(mr,4);

% variables for plotting
colIdx = [1 2 0 0 3 4];
vmr = reshape(vmr,[Ntrials*Nblock Nsubj]);
mr = reshape(mr,[Ntrials*Nblock Nsubj]);

% This section generates "alignment_matrix.csv" which contains the 
% relevant values of the alignment matrices for statistical analysis 
% in R.
% z = [vmr(1,:)'; vmr(40,:)'; vmr(41,:)'; mr(1,:)'; mr(40,:)'; mr(41,:)'];
% dlmwrite('alignment_matrix.csv',z);

figure(7); clf
% plot responses for all tracking blocks
for k = 1:2
    subplot(1,2,k); hold on 
    rectangle('Position',[9 -10 31 110],'FaceColor',[0 0 0 0.1],...
        'EdgeColor','none');
    plot([0 Ntrials*Nblock+1],[0 0],'k')
    
    % mean across subjects
    for i = 1:6
        idx = Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials;
        if i == 3 || i == 4 % plot blocks between Early and Late in black
            if k == 1
                plot(idx,vmr(idx,:),'k','Color',[0 0 0 0.5],'LineWidth',0.2)
                plot(idx,mean(vmr(idx,:),2),'k','LineWidth',1.5)
            else
                plot(idx,mr(idx,:),'k','Color',[0 0 0 0.5],'LineWidth',...
                    0.2)
                plot(idx,mean(mr(idx,:),2),'k','LineWidth',1.5)
            end
        else % all other blocks plot with colors
            if k == 1
                plot(idx,vmr(idx,:),'k','Color',[0 0 0 0.5],'LineWidth',...
                    0.2)
                plot(idx,mean(vmr(idx,:),2),'Color',col(colIdx(i),:)...
                    ,'LineWidth',1.5)
            else
                plot(idx,mr(idx,:),'k','Color',[0 0 0 0.5],'LineWidth',0.2)
                plot(idx,mean(mr(idx,:),2),'Color',col(colIdx(i),:)...
                    ,'LineWidth',1.5)
            end
        end
    end
    axis([0.5 Ntrials*Nblock+.5 -.2 0.8])
    set(gca,'TickDir','out')
    xticks(1:8:41)
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

% plot left panel of Figure 3C
figure(8); clf
subplot(1,2,1); hold on
rectangle('Position',[9 -10 31 110],'FaceColor',[0 0 0 0.1],'EdgeColor',...
    'none');
plot([0 Nblock*Ntrials+1],[0 0],'--k','LineWidth',1) % ideal baseline response
plot([0 Nblock*Ntrials+1],[90 90],'--k','LineWidth',1) % ideal compensation
% plot mean across subjects
for i = 1:Nblock
    plotIdx = Ntrials*(i-1)+1:Ntrials*(i-1)+Ntrials;
    plot(plotIdx,theta_opt2(plotIdx,:),'Color',[0 0 0 0.5],'LineWidth',0.2) % plot individual subjects
    if i == 3 || i == 4 % plot blocks between Early and Late in black
        plot(plotIdx,thetaMu(plotIdx),'Color',[0 0 0],'LineWidth',1.5)
    else % all other blocks plotted in color
        plot(plotIdx,thetaMu(plotIdx),'Color',colors(idx(i),:),...
            'LineWidth',1.5)
    end
end
set(gca,'TickDir','out')
title('Rotation')
xticks(1:8:41)
xlabel('Trial Number')
yticks(0:30:90)
ylabel(['Angle (' char(176) ')'])
axis([0.5 length(thetaMu) + 0.5 -10 100])

% plot right panel of Figure 3C
subplot(1,2,2); hold on
rectangle('Position',[9 -10 31 110],'FaceColor',[0 0 0 0.1],'EdgeColor',...
    'none');
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
        plot(plotIdx,scaleOrthMu(plotIdx),'Color',colors(idx(i),:),...
            'LineWidth',1.5)
    end
end
set(gca,'TickDir','out')
title('Mirror-Reversal')
xticks(1:8:41)
xlabel('Trial Number')
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