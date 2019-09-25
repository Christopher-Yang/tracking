clear all
load dat
groups = {'rot','rot_i'};
subj_rot = {'subj17','subj18','subj21','subj22','subj24','subj25','subj28','subj31','subj32','subj33'};
subj_rot_i = {'subj14','subj15','subj16','subj19','subj23','subj26','subj27','subj29','subj30','subj34'};
blocks = {'no_rot1','rot1','rot2','rot3','rot4','no_rot2'};
names = {'x_x_all','x_y_all','y_x_all','y_y_all'};
Nblock = length(blocks);
Nsubj = length(subj_rot);
Nfreq = length(data.rot.avg.x_x.freqs);
Ngroups = length(groups);
paramsInit = zeros([2*Nblock 1]);

for q = 1:2 % loop over groups
    if q == 1
        subj_name = subj_rot;
    else
        subj_name = subj_rot_i;
    end
    for p = 1:Nsubj % loop over subjects
        for k = 1:Nblock
            for i = 1:4
                cplx(:,i,k,p,q) = mean(data.(groups{q}).(subj_name{p}).(blocks{k}).(names{i}).ratio,2); % put all complex ratios in cplx
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

%% display results of fitting process
g = 2;
subj = 9;
blockIdx = 5;
freq = 1;

if g == 1
    subj_name = subj_rot;
else
    subj_name = subj_rot_i;
end

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
plot_subj(data.(groups{g}).(subj_name{subj}).(blocks{blockIdx}))
for i = 1:4
    subplot(2,2,gblocks(i)); hold on
    plot(phasor(i,:),'-ok','LineWidth',1.5,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none')
    plot([phasor(i,freq) cplx(freq,i,blockIdx,subj,g)],'b','LineWidth',2)
    if i <=2
        plot(dat([1 2],:),'.b','MarkerSize',10)
        if thetaOpt(i,blockIdx,freq) >= 0
            plot([0 real(slope(1))],[0 imag(slope(1))],'b','LineWidth',2)
        else
            plot([0 -real(slope(1))],[0 -imag(slope(1))],'b','LineWidth',2)
        end
    else
        plot(dat([3 4],:),'.b','MarkerSize',10)
        if thetaOpt(i,blockIdx,freq) >= 0
            plot([0 real(slope(2))],[0 imag(slope(2))],'b','LineWidth',2)
        else
            plot([0 -real(slope(2))],[0 -imag(slope(2))],'b','LineWidth',2)
        end
    end
    axis([-1.25 1.25 -1.25 1.25])
    pbaspect([1 1 1])
end

%% plot vectors and gain matrices
% for single subjects
% subj = 5;
% rotMat_mu = squeeze(rotMat(:,:,:,:,subj,:));

% for averaging across subjects
rotMat_mu = squeeze(mean(rotMat,5));

gblocks = [1 2 5 6];
figure(1); clf
for k = 1:4
    for i = 1:Nfreq
        subplot(4,7,7*(k-1)+i); hold on
        plot([0 rotMat_mu(1,1,i,gblocks(k),1)],[0 rotMat_mu(2,1,i,gblocks(k),1)],'LineWidth',1.5)
        plot([0 rotMat_mu(1,2,i,gblocks(k),1)],[0 rotMat_mu(2,2,i,gblocks(k),1)],'LineWidth',1.5)
        plot([0 1],[0 0],'k')
        plot([0 0],[0 1],'k')
        if i == 1
            ylabel('Baseline')
        end
        axis([-0.45 1 -0.45 1])
        axis square
    end
end

figure(2); clf
for k = 1:4
    for i = 1:Nfreq
        subplot(4,7,7*(k-1)+i); hold on
        plot([0 rotMat_mu(1,1,i,gblocks(k),2)],[0 rotMat_mu(2,1,i,gblocks(k),2)],'LineWidth',1.5)
        plot([0 rotMat_mu(1,2,i,gblocks(k),2)],[0 rotMat_mu(2,2,i,gblocks(k),2)],'LineWidth',1.5)
        plot([0 1],[0 0],'k')
        plot([0 0],[0 1],'k')
        if i == 1
            ylabel('Baseline')
        end
        axis([-0.45 1 -0.45 1])
        axis square
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gain matrices
col1 = [1 0 0];
col2 = [1 1 1];
Nstep = 100;
map1 = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2),Nstep)', linspace(col1(3),col2(3),Nstep)'];

col1 = [1 1 1];
col2 = [0 0 1];
map2 = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2),Nstep)', linspace(col1(3),col2(3),Nstep)'];

map = [map1; map2];
clims = [-1 1];

figure(3); clf
for k = 1:4
    for i = 1:Nfreq
        subplot(4,7,7*(k-1)+i)
        imagesc(rotMat_mu(:,:,i,gblocks(k),1),clims)
        colormap(map)
        set(gca,'TickDir','out','Xtick',[],'Ytick',[])
        axis square
        if i == 1
            ylabel('Baseline')
        end
    end
end

figure(4); clf
for k = 1:4
    for i = 1:Nfreq
        subplot(4,7,7*(k-1)+i)
        imagesc(rotMat_mu(:,:,i,gblocks(k),2),clims)
        colormap(map)
        set(gca,'TickDir','out','Xtick',[],'Ytick',[])
        axis square
        if i == 1
            ylabel('Baseline')
        end
    end
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
%     if i == 1
%         title
%     end
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