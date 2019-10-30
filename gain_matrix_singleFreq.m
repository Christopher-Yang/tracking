% clear all
% load dat
groups = {'rot','mir'};
blocks = {'no_rot1','rot1','rot2','rot3','rot4','no_rot2'};
graph_name = {'Baseline','Early','Train2','Train3','Late','Post'};
names = {'x_x_all','x_y_all','y_x_all','y_y_all'};
Nblock = length(blocks);
Nsubj = length(data.rot)-1;
Nfreq = length(data.rot{end}.freqX);
Ngroups = length(groups);
paramsInit = zeros([2*Nblock 1]);
labels = {'X_TX_H','X_TY_H','Y_TX_H','Y_TY_H'};

for q = 1:2 % loop over groups
    for p = 1:Nsubj % loop over subjects
        for k = 1:Nblock
            for i = 1:4
                cplx(:,i,k,p,q) = mean(data.(groups{q}){p}.(blocks{k}).phasors.Rhand.(names{i}).ratio,2); % put all complex ratios in cplx
            end
        end
        
        for i = 1:Nfreq
            x = cos(angle(cplx(i,[1 4],1,p,q)));
            y = sin(angle(cplx(i,[1 4],1,p,q)));
            template(:,i,p,q) = [(x(1)+y(1)*1j); (x(2)+y(2)*1j)]; % generate template of gain = 1 and phase same as baseline
            error = @(params) scale(params,cplx(i,[1 2],:,p,q),template(1,i,p,q));
            opt1(:,i,p,q) = fmincon(error,paramsInit); % fit gains for x frequencies
            error = @(params) scale(params,cplx(i,[3 4],:,p,q),template(2,i,p,q));
            opt2(:,i,p,q) = fmincon(error,paramsInit); % fit gains for y frequencies
        end
    end
end

% combine opt1 and opt2
thetaOpt = [reshape(opt1,[2 Nblock Nfreq Nsubj Ngroups]); reshape(opt2,[2 Nblock Nfreq Nsubj Ngroups])];

% shape thetaOpt into gain matrix format
rotMat = reshape(thetaOpt,[2 2 Nblock Nfreq Nsubj Ngroups]);
rotMat = permute(rotMat,[1 2 4 3 5 6]);

% fit theta to rotation matrices
for p = 1:Nsubj
    for k = 1:Nblock
        for i = 1:Nfreq
            thetaInit = 0;
            err2 = @(theta) fit_rotMat(theta,rotMat(:,:,i,k,p,1));
            thetaFit(i,k,p) = fmincon(err2,thetaInit);
        end
    end
end

% for plotting gain matrices
col1 = [1 0 0];
col2 = [1 1 1];
Nstep = 100;
map1 = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2),Nstep)', linspace(col1(3),col2(3),Nstep)'];

col1 = [1 1 1];
col2 = [0 0 1];
map2 = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2),Nstep)', linspace(col1(3),col2(3),Nstep)'];

map = [map1; map2];
clims = [-1 1];

% for plotting vectors
col1 = [0 128 0]/255;
col2 = [152 251 152]/255;
map1 = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2),Nfreq)', linspace(col1(3),col2(3),Nfreq)'];

col1 = [128 0 128]/255;
col2 = [230 230 250]/255;
map2 = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2),Nfreq)', linspace(col1(3),col2(3),Nfreq)'];

% mat2 = thetaOpt;
% mat2(3,:,:,:,1) = -mat2(3,:,:,:,1);
% mat2 = mean(mat2(2:3,[1 5 6],:,:,:),1);
% offAll = permute(mat2,[4 3 2 5 1]);
% z = reshape(offAll,[420 1]);
% dlmwrite('gain_matrix.csv',z);

%% display results of fitting process
g = 1;
subj = 4;
blockIdx = 5;
freq = 1;

% slope is the line which all data for a given frequency projects onto
slope = template(:,freq,subj,g);
for i = 1:numel(slope)
    slope(i) = slope(i)*(2/norm(slope(i)));
end

% reconstruct fitted phasor
template2 = template(:,:,subj,g);
template2 = [repmat(template2(1,:),[2 1]); repmat(template2(2,:),[2 1])];
phasor = squeeze(thetaOpt(:,blockIdx,:,subj,g)).*template2;
dat = squeeze(cplx(freq,:,:,subj,g)); % extract data to plot

gblocks = [1 3 2 4];
plot_subj(data.(groups{g}){subj}.(blocks{blockIdx}).phasors.Rhand)
for i = 1:4
    subplot(2,2,gblocks(i)); hold on
    plot(phasor(i,:),'-ok','LineWidth',1.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none')
%     plot([phasor(i,freq) cplx(freq,i,blockIdx,subj,g)],'b','LineWidth',2)
%     if i <=2
%         plot(dat([1 2],:),'.b','MarkerSize',10)
%         if thetaOpt(i,blockIdx,freq) >= 0
%             plot([0 real(slope(1))],[0 imag(slope(1))],'b','LineWidth',2)
%         else
%             plot([0 -real(slope(1))],[0 -imag(slope(1))],'b','LineWidth',2)
%         end
%     else
%         plot(dat([3 4],:),'.b','MarkerSize',10)
%         if thetaOpt(i,blockIdx,freq) >= 0
%             plot([0 real(slope(2))],[0 imag(slope(2))],'b','LineWidth',2)
%         else
%             plot([0 -real(slope(2))],[0 -imag(slope(2))],'b','LineWidth',2)
%         end
%     end
    axis([-1.25 1.25 -1.25 1.25])
    pbaspect([1 1 1])
end

%% plot vectors and gain matrices
% for single subjects
% subj = 6;
% mat = squeeze(thetaOpt(:,:,:,subj,:));

% for averaging across subjects
mat = squeeze(mean(thetaOpt,4));

gblocks = [1 2 5 6];
figure(1); clf
for i = 1:length(gblocks)
    subplot(2,length(gblocks),i); hold on
    plot([0 1],[0 0],'k')
    plot([0 0],[0 1],'k')
    for k = 1:Nfreq
        plot([0 mat(1,gblocks(i),k,1)],[0 mat(2,gblocks(i),k,1)],'LineWidth',1,'Color',map1(k,:))
        plot([0 mat(3,gblocks(i),k,1)],[0 mat(4,gblocks(i),k,1)],'LineWidth',1,'Color',map2(k,:))
    end
    axis([-0.75 1 -0.75 1])
    axis square
    title(graph_name(gblocks(i)))
    if i == 1
        ylabel('Rotation')
    end
    
    subplot(2,length(gblocks),i+length(gblocks)); hold on
    plot([0 1],[0 0],'k')
    plot([0 0],[0 1],'k')
    for k = 1:Nfreq
        plot([0 mat(1,gblocks(i),k,2)],[0 mat(2,gblocks(i),k,2)],'LineWidth',1,'Color',map1(k,:))
        plot([0 mat(3,gblocks(i),k,2)],[0 mat(4,gblocks(i),k,2)],'LineWidth',1,'Color',map2(k,:))
    end
    axis([-0.75 1 -0.75 1])
    axis square
    if i == 1
        ylabel('Mirror-Reversal')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gain matrices
gblocks = 1:6;
figure(2); clf
for i = 1:Nfreq
    subplot(2,Nfreq,i)
    imagesc(mat(:,gblocks,i,1),clims)
    colormap(map)
    set(gca,'TickDir','out','Xtick',[],'Ytick',[])
    if i == 1
        title('Low freq')
        ylabel('Rotation')
        yticks(1:4)
        yticklabels(labels)
    elseif i == Nfreq
        title('High freq')
    end

    subplot(2,Nfreq,i+Nfreq)
    imagesc(mat(:,gblocks,i,2),clims)
    colormap(map)
    set(gca,'TickDir','out','Xtick',[],'Ytick',[])
    if i == 1
        title('Low freq')
        ylabel('Mirror-Reversal')
        yticks(1:4)
        yticklabels(labels)
    elseif i == Nfreq
        title('High freq')
    end
end

%% plot mirror-reversal compensation perpendicular to mirroring axis
gblocks = [1:2 5:6];
colors = lines;
colors = colors(1:7,:);

R = rotz(-45);
R = R(1:2,1:2);
perpAxis = R*[1 0]';

for k = 1:Nsubj
    for p = 1:Nblock
        for q = 1:Nfreq
            comp(q,p,k) = perpAxis'*rotMat(:,:,q,p,k,2)*perpAxis;
        end
    end
end

comp2 = mean(comp,3);
compSE = std(comp,[],3)/sqrt(Nsubj);
thetaFitMu = mean(thetaFit,3);
thetaFitSE = std(thetaFit,[],3)/sqrt(Nsubj);

figure(3); clf
subplot(1,2,1); hold on
for i = 1:length(gblocks)
    errorbar(thetaFitMu(:,gblocks(i)),thetaFitSE(:,gblocks(i)),'-o','MarkerFaceColor',colors(i,:),'MarkerEdgeColor','none','LineWidth',1)
end
plot([1 Nfreq],[0 0],'k','LineWidth',1)
plot([1 Nfreq],[90 90],'--k','LineWidth',1)
title('Rotation')
xticks([1 Nfreq])
xticklabels({'Low Freq','High Freq'})
yticks(0:30:90)
ylabel('Fitted angle')
axis([1 Nfreq -10 90])

subplot(1,2,2); hold on
for i = 1:length(gblocks)
    errorbar(comp2(:,gblocks(i)),compSE(:,gblocks(i)),'-o','MarkerFaceColor',colors(i,:),'MarkerEdgeColor','none','LineWidth',1)
end
plot([1 Nfreq],[1 1],'k','LineWidth',1)
plot([1 Nfreq],[-1 -1],'--k','LineWidth',1)
title('Mirror-Reversal')
xticks([1 Nfreq])
xticklabels({'Low Freq','High Freq'})
yticks(-1:0.5:1)
ylabel('Gain orthogonal to mirroring axis')
axis([1 Nfreq -1 1])
legend(graph_name(gblocks),'Location','southeast')

%%
col = copper;
colors = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

mat2 = thetaOpt;
mat2(3,:,:,:,1) = -mat2(3,:,:,:,1);
mat2 = [mean(mat2([1 4],:,:,:,:),1); mean(mat2(2:3,:,:,:,:),1)];
on = squeeze(mean(mat2(1,:,:,:,:),4));
off = squeeze(mean(mat2(2,:,:,:,:),4));

onAll = permute(mat2(1,:,:,:,:),[2 4 3 5 1]);
offAll = permute(mat2(2,:,:,:,:),[2 4 3 5 1]);

onSE = squeeze(std(onAll,[],2))/sqrt(Nsubj);
offSE = squeeze(std(offAll,[],2))/sqrt(Nsubj);

gblocks = 1:6;
figure(4); clf
for i = 1:Nfreq
    subplot(2,2,1); hold on
    plot([1 length(gblocks)],[0 0],'--k','LineWidth',1)
    plot(on(gblocks,i,1),'-o','Color',colors(i,:),'LineWidth',1,'MarkerFaceColor',colors(i,:),'MarkerEdgeColor','none')
    errorbar(on(gblocks,i,1),onSE(gblocks,i,1),'-o','Color',colors(i,:),'LineWidth',1,'MarkerFaceColor',colors(i,:),'MarkerEdgeColor','none')
    title('Rotation')
    ylim([-0.25 0.8])
    xticks(1:length(gblocks))
    xticklabels(graph_name(gblocks))
    if i == 1
        ylabel('On-axis')
    end
    pbaspect([1 1 1])
    
    subplot(2,2,3); hold on
    plot([1 length(gblocks)],[0 0],'--k','LineWidth',1)
%     plot(off(gblocks,i,1),'-o','Color',colors(i,:),'LineWidth',1,'MarkerFaceColor',colors(i,:),'MarkerEdgeColor','none')
    errorbar(off(gblocks,i,1),offSE(gblocks,i,1),'-o','Color',colors(i,:),'LineWidth',1,'MarkerFaceColor',colors(i,:),'MarkerEdgeColor','none')
    ylim([-0.1 0.75])
    yticks(0:0.3:0.6)
    xticks(1:length(gblocks))
    xticklabels(graph_name(gblocks))
    if i == 1
        ylabel('Cross-axis gain')
    end
    pbaspect([1 1 1])

    subplot(2,2,2); hold on
    plot([1 length(gblocks)],[0 0],'--k','LineWidth',1)
%     plot(on(gblocks,i,2),'-o','Color',colors(i,:),'LineWidth',1,'MarkerFaceColor',colors(i,:),'MarkerEdgeColor','none')
    errorbar(on(gblocks,i,2),onSE(gblocks,i,2),'-o','Color',colors(i,:),'LineWidth',1,'MarkerFaceColor',colors(i,:),'MarkerEdgeColor','none')
    title('Mirror-Reversal')
    ylim([-0.25 0.8])
    xticks(1:length(gblocks))
    xticklabels(graph_name(gblocks))
    pbaspect([1 1 1])
    
    subplot(2,2,4); hold on
    plot([1 length(gblocks)],[0 0],'--k','LineWidth',1)
%     plot(off(gblocks,i,2),'-o','Color',colors(i,:),'LineWidth',1,'MarkerFaceColor',colors(i,:),'MarkerEdgeColor','none')
    errorbar(off(gblocks,i,2),offSE(gblocks,i,2),'-o','Color',colors(i,:),'LineWidth',1,'MarkerFaceColor',colors(i,:),'MarkerEdgeColor','none')
    ylim([-0.1 0.75])
    yticks(0:0.3:0.6)
    xticks(1:length(gblocks))
    xticklabels(graph_name(gblocks))
    pbaspect([1 1 1])
end

%% graph average fitted phasors across subjects
col1 = [1 0.9 0.3];
col2 = [1 0 0];
colors = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2),Nfreq)', linspace(col1(3),col2(3),Nfreq)'];
gblocks = [1 3 2 4];

group = 1;
blockIdx = 6;

phasorMat = repmat(permute(phasor(:,:,group),[1 3 2]),[1 4]);
phasorLines = permute(thetaOpt(:,blockIdx,:,:,group),[3 1 4 2]).*phasorMat;
phasorLines = mean(phasorLines,3);

dat_all = permute(squeeze(cplx(:,:,blockIdx,:,group)),[3 1 2]);
dat = squeeze(mean(dat_all,1));

figure(6); clf
for k = 1:4
    subplot(2,2,gblocks(k)); hold on
    for i = 1:Nfreq
        plot(dat_all(:,i,k),'.','Color',colors(i,:),'MarkerSize',10)
        if i > 1
            plot(dat(i-1:i,k),'-o','Color',colors(i,:),'LineWidth',1,'MarkerFaceColor',colors(i-1,:),'MarkerEdgeColor','none')
        end
    end
    plot(phasorLines(:,k),'-ok','LineWidth',2,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none')
    plot([-1.5 1.5],[0 0],'k','LineWidth',1)
    plot([0 0],[-1.5 1.5],'k','LineWidth',1)
    axis([-1.25 1.25 -1.25 1.25])
    pbaspect([1 1 1])
end

%%
function e = scale(params,cplx_data,template)
    theta = reshape(params,[2 length(params)/2]);
    phasors = theta*template;
    e = (phasors - squeeze(cplx_data)).^2;
    for k = 1:numel(e)
        e(k) = norm(e(k));
    end
    e = sum(sum(e));
end

function e = fit_rotMat(theta,rotMat_opt)
    rot = rotz(theta);
    rot = rot(1:2,1:2);
    e = (rot-rotMat_opt).^2;
    e = sum(sum(e));
end