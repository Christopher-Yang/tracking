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
        
        z = zeros([15 1]);
        
        for i = 1:7
            paramsInit = [cplx(i,1,1); z];
            error = @(params) scale(params,cplx(i,:,:));
            paramsOpt = fmincon(error,paramsInit);
            phasor(i,p,q) = paramsOpt(1);
            thetaOpt = [1; paramsOpt(2:end)];
            thetaOpt = reshape(thetaOpt, [4 4]);
        
            for k = 1:length(blocks)
                rotMat(:,:,i,k,p,q) = reshape(thetaOpt(:,k), [2 2]);
                if q == 1
                    thetaInit = 0;
                    err2 = @(theta) fit_rotMat(theta,rotMat(:,:,i,k,p,1));
                    thetaFit(i,k,p) = fmincon(err2,thetaInit);
                end
            end
            
            ph = abs(phasor(i,p,q));
            lambda = mean(ph);
        
            if q == 1
                vmrBase(:,:,i,p) = rotMat(:,:,i,1,p,q).*lambda;
                vmrEarly(:,:,i,p) = rotMat(:,:,i,2,p,q).*lambda;
                vmrLate(:,:,i,p) = rotMat(:,:,i,3,p,q).*lambda;
                vmrAfter(:,:,i,p) = rotMat(:,:,i,4,p,q).*lambda;
            else
                mrBase(:,:,i,p) = rotMat(:,:,i,1,p,q).*lambda;
                mrEarly(:,:,i,p) = rotMat(:,:,i,2,p,q).*lambda;
                mrLate(:,:,i,p) = rotMat(:,:,i,3,p,q).*lambda;
                mrAfter(:,:,i,p) = rotMat(:,:,i,4,p,q).*lambda;
            end
        end
    end
end

thetaFit_bar = mean(thetaFit,3);
thetaFit_se = std(thetaFit,[],3)/sqrt(10);

%% plot vectors
% for single subjects
% subj = 5;
% vmrBaseMu = vmrBase(:,:,:,subj);
% vmrEarlyMu = vmrEarly(:,:,:,subj);
% vmrLateMu = vmrLate(:,:,:,subj);
% vmrAfterMu = vmrAfter(:,:,:,subj);
% 
% mrBaseMu = mrBase(:,:,:,subj);
% mrEarlyMu = mrEarly(:,:,:,subj);
% mrLateMu = mrLate(:,:,:,subj);
% mrAfterMu = mrAfter(:,:,:,subj);

% for averaging across subjects
vmrBaseMu = mean(vmrBase,4);
vmrEarlyMu = mean(vmrEarly,4);
vmrLateMu = mean(vmrLate,4);
vmrAfterMu = mean(vmrAfter,4);

mrBaseMu = mean(mrBase,4);
mrEarlyMu = mean(mrEarly,4);
mrLateMu = mean(mrLate,4);
mrAfterMu = mean(mrAfter,4);

figure(1); clf
for i = 1:7
    subplot(4,7,i); hold on
    plot([0 vmrBaseMu(1,1,i)],[0 vmrBaseMu(2,1,i)],'LineWidth',1.5)
    plot([0 vmrBaseMu(1,2,i)],[0 vmrBaseMu(2,2,i)],'LineWidth',1.5)
    plot([0 1],[0 0],'k')
    plot([0 0],[0 1],'k')
    if i == 1
        ylabel('Baseline')
    end
    axis([-0.45 1 -0.45 1])
    axis square
    
    subplot(4,7,i+7); hold on
    plot([0 vmrEarlyMu(1,1,i)],[0 vmrEarlyMu(2,1,i)],'LineWidth',1.5)
    plot([0 vmrEarlyMu(1,2,i)],[0 vmrEarlyMu(2,2,i)],'LineWidth',1.5)
    plot([0 1],[0 0],'k')
    plot([0 0],[0 1],'k')
    if i == 1
        ylabel('Early')
    end
    axis([-0.45 1 -0.45 1])
    axis square
    
    subplot(4,7,i+14); hold on
    plot([0 vmrLateMu(1,1,i)],[0 vmrLateMu(2,1,i)],'LineWidth',1.5)
    plot([0 vmrLateMu(1,2,i)],[0 vmrLateMu(2,2,i)],'LineWidth',1.5)
    plot([0 1],[0 0],'k')
    plot([0 0],[0 1],'k')
    if i == 1
        ylabel('Late')
    end
    axis([-0.45 1 -0.45 1])
    axis square
    
    subplot(4,7,i+21); hold on
    plot([0 vmrAfterMu(1,1,i)],[0 vmrAfterMu(2,1,i)],'LineWidth',1.5)
    plot([0 vmrAfterMu(1,2,i)],[0 vmrAfterMu(2,2,i)],'LineWidth',1.5)
    plot([0 1],[0 0],'k')
    plot([0 0],[0 1],'k')
    if i == 1
        ylabel('Post')
    end
    axis([-0.45 1 -0.45 1])
    axis square
end

figure(2); clf
for i = 1:7
    subplot(4,7,i); hold on
    plot([0 mrBaseMu(1,1,i)],[0 mrBaseMu(2,1,i)],'LineWidth',1.5)
    plot([0 mrBaseMu(1,2,i)],[0 mrBaseMu(2,2,i)],'LineWidth',1.5)
    plot([0 1],[0 0],'k')
    plot([0 0],[0 1],'k')
    if i == 1
        ylabel('Baseline')
    end
    axis([-0.45 1 -0.45 1])
    axis square
    
    subplot(4,7,i+7); hold on
    plot([0 mrEarlyMu(1,1,i)],[0 mrEarlyMu(2,1,i)],'LineWidth',1.5)
    plot([0 mrEarlyMu(1,2,i)],[0 mrEarlyMu(2,2,i)],'LineWidth',1.5)
    plot([0 1],[0 0],'k')
    plot([0 0],[0 1],'k')
    if i == 1
        ylabel('Early')
    end
    axis([-0.45 1 -0.45 1])
    axis square
    
    subplot(4,7,i+14); hold on
    plot([0 mrLateMu(1,1,i)],[0 mrLateMu(2,1,i)],'LineWidth',1.5)
    plot([0 mrLateMu(1,2,i)],[0 mrLateMu(2,2,i)],'LineWidth',1.5)
    plot([0 1],[0 0],'k')
    plot([0 0],[0 1],'k')
    if i == 1
        ylabel('Late')
    end
    axis([-0.45 1 -0.45 1])
    axis square
    
    subplot(4,7,i+21); hold on
    plot([0 mrAfterMu(1,1,i)],[0 mrAfterMu(2,1,i)],'LineWidth',1.5)
    plot([0 mrAfterMu(1,2,i)],[0 mrAfterMu(2,2,i)],'LineWidth',1.5)
    plot([0 1],[0 0],'k')
    plot([0 0],[0 1],'k')
    if i == 1
        ylabel('Late')
    end
    axis([-0.45 1 -0.45 1])
    axis square
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
for i = 1:7
    subplot(4,7,i)
    imagesc(vmrBaseMu(:,:,i),clims)
    colormap(map)
    set(gca,'TickDir','out','Xtick',[],'Ytick',[])
    axis square
    if i == 1
        ylabel('Baseline')
    end

    subplot(4,7,i+7)
    imagesc(vmrEarlyMu(:,:,i),clims)
    colormap(map)
    set(gca,'TickDir','out','Xtick',[],'Ytick',[])
    axis square
    if i == 1
        ylabel('Early')
    end

    subplot(4,7,i+14)
    imagesc(vmrLateMu(:,:,i),clims)
    colormap(map)
    set(gca,'TickDir','out','Xtick',[],'Ytick',[])
    axis square
    if i == 1
        ylabel('Late')
    end
    
    subplot(4,7,i+21)
    imagesc(vmrAfterMu(:,:,i),clims)
    colormap(map)
    set(gca,'TickDir','out','Xtick',[],'Ytick',[])
    axis square
    if i == 1
        ylabel('After')
    end
end

figure(4); clf
for i = 1:7
    subplot(4,7,i)
    imagesc(mrBaseMu(:,:,i),clims)
    colormap(map)
    set(gca,'TickDir','out','Xtick',[],'Ytick',[])
    axis square
    if i == 1
        ylabel('Baseline')
    end
    
    subplot(4,7,i+7)
    imagesc(mrEarlyMu(:,:,i),clims)
    colormap(map)
    set(gca,'TickDir','out','Xtick',[],'Ytick',[])
    if i == 1
        ylabel('Early')
    end
    axis square
    
    subplot(4,7,i+14)
    imagesc(mrLateMu(:,:,i),clims)
    colormap(map)
    set(gca,'Xtick',[],'Ytick',[])
    if i == 1
        ylabel('Late')
    end
    axis square
    
    subplot(4,7,i+21)
    imagesc(mrAfterMu(:,:,i),clims)
    colormap(map)
    set(gca,'Xtick',[],'Ytick',[])
    if i == 1
        ylabel('Post')
    end
    axis square
end

%% fitted phasors for individual subjects
group = 1;
subj = 3;
blockIdx = 4;

for i = 1:7
    line11(i) = rotMat(1,1,i,blockIdx,subj,group)*phasor(i,subj,group);
    line12(i) = rotMat(1,2,i,blockIdx,subj,group)*phasor(i,subj,group);
    line21(i) = rotMat(2,1,i,blockIdx,subj,group)*phasor(i,subj,group);
    line22(i) = rotMat(2,2,i,blockIdx,subj,group)*phasor(i,subj,group);
end

plot_subj(data.rot.(subj_rot{subj}).(blocks{blockIdx}),5)
for i = 1:7
    subplot(2,2,1)
    plot(line11,'k','LineWidth',3)
    subplot(2,2,2)
    plot(line12,'k','LineWidth',3)
    subplot(2,2,3)
    plot(line21,'k','LineWidth',3)
    subplot(2,2,4)
    plot(line22,'k','LineWidth',3)
end

%%
function e = scale(params,cplx_data)
    cplx_fit = params(1);
    theta = [1; params(2:end)];
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

function e = fit_rotMat(theta,rotMat_opt)
    rot = rotz(theta);
    rot = rot(1:2,1:2);
    e = (rot-rotMat_opt).^2;
    e = sum(sum(e));
end