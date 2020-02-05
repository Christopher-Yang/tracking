function graph_gainMatrix(data,experiment)

% set variables for analysis
groups = {'rot','mir'};
blocks = {'baseline','pert1','pert2','pert3','pert4','post'};
graph_name = {'Baseline','Early','Train2','Train3','Late','Post'};
names = {'x_x_all','x_y_all','y_x_all','y_y_all'};
Nblock = length(blocks);
Nsubj = length(data.rot)-1;
Nfreq = length(data.rot{end}.freqX);
Ngroup = length(groups);

paramsInit = zeros([2*Nblock 1]); % initialize parameters

% perform fitting
for q = 1:Ngroup % loop over groups
    for p = 1:Nsubj % loop over subjects
        for k = 1:Nblock % loop over blocks
            for i = 1:4 % loop over combinations of target and hand movement axes
                % i=1: x-target to x-hand
                % i=2: x-target to y-hand
                % i=3: y-target to x-hand
                % i=4: y-target to y-hand
                cplx(:,i,k,p,q) = mean(data.(groups{q}){p}.(blocks{k}).phasors.Rhand.(names{i}).ratio,2); % put all complex ratios in cplx
            end
        end
        
        for i = 1:Nfreq
            x = cos(angle(cplx(i,[1 4],1,p,q))); % x-component of xx and yy baseline phasors
            y = sin(angle(cplx(i,[1 4],1,p,q))); % y-component of xx and yy baseline phasors
            template(:,i,p,q) = [(x(1)+y(1)*1j); (x(2)+y(2)*1j)]; % generate template of gain = 1 and phase same as baseline
            error = @(params) scale(params,cplx(i,[1 2],:,p,q),template(1,i,p,q));
            opt1(:,i,p,q) = fmincon(error,paramsInit); % fit gains for x frequencies
            error = @(params) scale(params,cplx(i,[3 4],:,p,q),template(2,i,p,q));
            opt2(:,i,p,q) = fmincon(error,paramsInit); % fit gains for y frequencies
        end
    end
end

thetaOpt = [reshape(opt1,[2 Nblock Nfreq Nsubj Ngroup]); reshape(opt2,[2 Nblock Nfreq Nsubj Ngroup])]; % combine opt1 and opt2
gainMat = reshape(thetaOpt,[2 2 Nblock Nfreq Nsubj Ngroup]); % shape thetaOpt into gain matrix format
gainMat = permute(gainMat,[1 2 4 3 5 6]);

% fit theta to rotation matrices
for p = 1:Nsubj
    for k = 1:Nblock
        for i = 1:Nfreq
            thetaInit = 0;
            err2 = @(theta) fit_rotMat(theta,gainMat(:,:,i,k,p,1));
            thetaFit(i,k,p) = fmincon(err2,thetaInit);
            while thetaFit(i,k,p) >= 180
                thetaFit(i,k,p) = thetaFit(i,k,p) - 360;
            end
            while thetaFit(i,k,p) < -180
                thetaFit(i,k,p) = thetaFit(i,k,p) + 360;
            end 
        end
    end
end

% compute orthogonal gain for mirror-reversal group
R = rotz(-45); % 45 degree clockwise rotation matrix
R = R(1:2,1:2);
perpAxis = R*[1 0]';

for k = 1:Nsubj
    for p = 1:Nblock
        for q = 1:Nfreq
            gainOrth(q,p,k) = perpAxis'*gainMat(:,:,q,p,k,2)*perpAxis; 
        end
    end
end

%% plot vectors from gain matrices (Figure 5A, S2B, and S5B)
% generate color maps
col1 = [0 128 0]/255;
col2 = [152 251 152]/255;
map1 = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2),Nfreq)', linspace(col1(3),col2(3),Nfreq)'];

col1 = [128 0 128]/255;
col2 = [230 230 250]/255;
map2 = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2),Nfreq)', linspace(col1(3),col2(3),Nfreq)'];

gblocks = [1 2 5 6]; % blocks to plot

if experiment == 1 % figures for main experiment (5A and S2B)
% This loop makes two figures. p=1 generates Figure 5A for and p=2 generates
% Figure S2B.
    for p = 1:2
        figure(p+15); clf
        if p == 1
            mat = squeeze(mean(thetaOpt,4)); % make plot for average of both rotation and mirror-reversal groups
        end
        for i = 1:length(gblocks) % iterate of the blocks to plot
            if p == 2
                mat = squeeze(thetaOpt(:,:,:,4,:)); % make plot for rotation group subj4
            end
            subplot(2,length(gblocks),i); hold on
            plot([0 1],[0 0],'k') % unit x vector
            plot([0 0],[0 1],'k') % unit y vector
            for k = 1:Nfreq
                plot([0 mat(1,gblocks(i),k,1)],[0 mat(2,gblocks(i),k,1)],'LineWidth',1,'Color',map1(k,:)) % plot green vectors
                plot([0 mat(3,gblocks(i),k,1)],[0 mat(4,gblocks(i),k,1)],'LineWidth',1,'Color',map2(k,:)) % plot purple vectors
            end
            axis([-0.75 1 -0.75 1])
            axis square
            title(graph_name(gblocks(i)))
            if i == 1
                ylabel('Rotation')
            end

            if p == 2
                mat = squeeze(thetaOpt(:,:,:,9,:)); % make plot for mirror-reversal group subj9
            end
            subplot(2,length(gblocks),i+length(gblocks)); hold on
            plot([0 1],[0 0],'k') % unit x vector
            plot([0 0],[0 1],'k') % unit y vector
            for k = 1:Nfreq
                plot([0 mat(1,gblocks(i),k,2)],[0 mat(2,gblocks(i),k,2)],'LineWidth',1,'Color',map1(k,:)) % plot green vectors
                plot([0 mat(3,gblocks(i),k,2)],[0 mat(4,gblocks(i),k,2)],'LineWidth',1,'Color',map2(k,:)) % plot purple vectors
            end
            axis([-0.75 1 -0.75 1])
            axis square
            if i == 1
                ylabel('Mirror-Reversal')
            end
        end
    end
else % figures for second experiment (S5B)
    figure(26); clf 
    mat = squeeze(mean(thetaOpt,4)); % make plot for average of both rotation and mirror-reversal groups
    for i = 1:length(gblocks) % iterate of the blocks to plot
        subplot(2,length(gblocks),i); hold on
        plot([0 1],[0 0],'k') % unit x vector
        plot([0 0],[0 1],'k') % unit y vector
        for k = 1:Nfreq
            plot([0 mat(1,gblocks(i),k,1)],[0 mat(2,gblocks(i),k,1)],'LineWidth',1,'Color',map1(k,:)) % plot green vectors
            plot([0 mat(3,gblocks(i),k,1)],[0 mat(4,gblocks(i),k,1)],'LineWidth',1,'Color',map2(k,:)) % plot purple vectors
        end
        axis([-0.75 1 -0.75 1])
        axis square
        title(graph_name(gblocks(i)))
        if i == 1
            ylabel('Rotation')
        end
        
        subplot(2,length(gblocks),i+length(gblocks)); hold on
        plot([0 1],[0 0],'k') % unit x vector
        plot([0 0],[0 1],'k') % unit y vector
        for k = 1:Nfreq
            plot([0 mat(1,gblocks(i),k,2)],[0 mat(2,gblocks(i),k,2)],'LineWidth',1,'Color',map1(k,:)) % plot green vectors
            plot([0 mat(3,gblocks(i),k,2)],[0 mat(4,gblocks(i),k,2)],'LineWidth',1,'Color',map2(k,:)) % plot purple vectors
        end
        axis([-0.75 1 -0.75 1])
        axis square
        if i == 1
            ylabel('Mirror-Reversal')
        end
    end
end

%% plot off-diagonal elements (Figure 5B and S5C)
col = copper;
colors = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

% average the off-diagonal elements together
mat2 = thetaOpt; 
mat2(3,:,:,:,1) = -mat2(3,:,:,:,1); % for the rotation group, flip the sign of element 1,2
mat2 = [mean(mat2([1 4],:,:,:,:),1); mean(mat2(2:3,:,:,:,:),1)]; % average on-diagonal and off-diagonal elements
on = squeeze(mean(mat2(1,:,:,:,:),4)); % extract on-diagonal mean
off = squeeze(mean(mat2(2,:,:,:,:),4)); % extract off-diagonal mean

% compute standard error over the averaged on- and off-diagonal elements
onAll = permute(mat2(1,:,:,:,:),[2 4 3 5 1]);
offAll = permute(mat2(2,:,:,:,:),[2 4 3 5 1]);
onSE = squeeze(std(onAll,[],2))/sqrt(Nsubj);
offSE = squeeze(std(offAll,[],2))/sqrt(Nsubj);

gblocks = 1:6; % blocks to plot
if experiment == 1 % figure for main experiment (5B)
    f = 18;
else % figure for second experiment (S5C)
    f = 27;
end

figure(f); clf
for i = 1:Nfreq
    subplot(1,2,1); hold on % rotation group
    plot([1 length(gblocks)],[0 0],'--k','LineWidth',1)
    errorbar(off(gblocks,i,1),offSE(gblocks,i,1),'-o','Color',colors(i,:),'LineWidth',1,'MarkerFaceColor',colors(i,:),'MarkerEdgeColor','none')
    axis([1 6 -0.1 0.75])
    yticks(0:0.3:0.6)
    xticks([1 2 5 6])
    xticklabels(graph_name([1 2 5 6]))
    title('Rotation')
    if i == 1
        ylabel('Cross-axis gain')
    end
    pbaspect([1 1 1])
    
    subplot(1,2,2); hold on % mirror-reversal group
    plot([1 length(gblocks)],[0 0],'--k','LineWidth',1)
    errorbar(off(gblocks,i,2),offSE(gblocks,i,2),'-o','Color',colors(i,:),'LineWidth',1,'MarkerFaceColor',colors(i,:),'MarkerEdgeColor','none')
    axis([1 6 -0.1 0.75])
    yticks(0:0.3:0.6)
    xticks([1 2 5 6])
    xticklabels(graph_name([1 2 5 6]))
    title('Mirror Reversal')
    pbaspect([1 1 1])
end

%% plot rotation angle and orthogonal gain (Figure 5C and S5D)
gblocks = [1:2 5:6];
colors = lines;
colors = colors(1:7,:);

gainOrthMu = mean(gainOrth,3);
gainOrthSE = std(gainOrth,[],3)/sqrt(Nsubj);
thetaFitMu = mean(thetaFit,3);
thetaFitSE = std(thetaFit,[],3)/sqrt(Nsubj);

if experiment == 1
    f = 19;
else
    f = 28;
end

figure(f); clf
subplot(1,2,1); hold on
plot([1 Nfreq],[0 0],'k','LineWidth',1) % ideal baseline response
plot([1 Nfreq],[90 90],'--k','LineWidth',1) % ideal compensation
for i = 1:length(gblocks)
    errorbar(thetaFitMu(:,gblocks(i)),thetaFitSE(:,gblocks(i)),'-o','MarkerFaceColor',colors(i,:),'MarkerEdgeColor','none','LineWidth',1)
end
title('Rotation')
xticks([1 Nfreq])
xticklabels({'Low Freq','High Freq'})
yticks(0:30:90)
ylabel('Fitted angle')
axis([1 Nfreq -10 90])

subplot(1,2,2); hold on
plot([1 Nfreq],[1 1],'k','LineWidth',1,'HandleVisibility','off') % ideal baseline response
plot([1 Nfreq],[-1 -1],'--k','LineWidth',1,'HandleVisibility','off') % ideal compensation
plot([1 Nfreq],[0 0],'k','HandleVisibility','off')
for i = 1:length(gblocks)
    errorbar(gainOrthMu(:,gblocks(i)),gainOrthSE(:,gblocks(i)),'-o','MarkerFaceColor',colors(i,:),'MarkerEdgeColor','none','LineWidth',1)
end
title('Mirror-Reversal')
xticks([1 Nfreq])
xticklabels({'Low Freq','High Freq'})
yticks(-1:0.5:1)
ylabel('Gain orthogonal to mirroring axis')
axis([1 Nfreq -1 1])
legend(graph_name(gblocks),'Location','southeast')

%% plot gain matrices (Figure S3)
if experiment == 1
    % generate color map
    col1 = [1 0 0];
    col2 = [1 1 1];
    Nstep = 100;
    map1 = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2),Nstep)', linspace(col1(3),col2(3),Nstep)'];
    
    col1 = [1 1 1];
    col2 = [0 0 1];
    map2 = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2),Nstep)', linspace(col1(3),col2(3),Nstep)'];
    
    map = [map1; map2];
    clims = [-1 1];
    
    % choose whether to make plot for single subject...
    % subj = 6;
    % mat = squeeze(thetaOpt(:,:,:,subj,:));
    
    % OR average across subjects
    mat = squeeze(mean(gainMat,5));
    
    % i=1 plots rotation data, i=2 plots mirror-reversal data
    for i = 1:Ngroup
        figure(i+19); clf
        for k = 1:Nblock
            for p = 1:Nfreq
                subplot(Nblock,Nfreq,(k-1)*Nfreq+p)
                imagesc(mat(:,:,p,k,i),clims)
                colormap(map)
                set(gca,'Xtick',[],'Ytick',[])
                pbaspect([1 1 1])
                if p == 1
                    switch k
                        case 1
                            ylabel('Baseline')
                        case 2
                            ylabel('Early')
                        case 5
                            ylabel('Late')
                        case 6
                            ylabel('Post')
                    end
                elseif p == 4
                    if k == Nblock
                        xlabel('Frequency')
                    elseif k == 1
                        if i == 1
                            title('Rotation')
                        else
                            title('Mirror Reversal')
                        end
                    end
                end
            end
        end
    end
end

%%
function e = scale(params,cplx_data,template)
    theta = reshape(params,[2 length(params)/2]);
    phasors = theta*template;
    e = (phasors - squeeze(cplx_data)).^2;
    for v = 1:numel(e)
        e(v) = norm(e(v));
    end
    e = sum(sum(e));
end

function e = fit_rotMat(theta,rotMat_opt)
    rot = rotz(theta);
    rot = rot(1:2,1:2);
    e = (rot-rotMat_opt).^2;
    e = sum(sum(e));
end
end