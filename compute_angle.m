Nsubj = length(data)-1;
Nfreq = 6;
names = {'x_x_all','x_y_all','y_x_all','y_y_all'};

for p = 1:Nsubj
    for k = 1:length(block_name)
        for i = 1:length(names)
            cplx(:,i,k) = mean(data{p}.(block_name{k}).phasors.cursor.(names{i}).ratio,2);
        end
    end
    
    theta = zeros([length(block_name)*4-1 1]);
    paramsInit = [cplx(:,1,1); theta];
    
    error = @(params) scale(params,cplx,Nfreq);
    paramsOpt = fmincon(error,paramsInit);
    phasor(:,p) = paramsOpt(1:Nfreq);
    thetaOpt = [1; paramsOpt(Nfreq+1:end)];
    Nblocks = length(thetaOpt)/4;
    thetaOpt = reshape(thetaOpt, [4 Nblocks]);
    ph = abs(phasor(:,p));
    lambda = mean(ph);
    
    for k = 1:length(thetaOpt)
%         rotMat(:,:,k,p) = reshape(thetaOpt(:,k), [2 2]);
        rotMat.(block_name{k})(:,:,p) = reshape(thetaOpt(:,k), [2 2]);
        rotMat.(block_name{k})(:,:,p) = rotMat.(block_name{k})(:,:,p).*lambda;
    end
    
%     Base(:,:,p) = rotMat(:,:,1,p).*lambda;
%     Early(:,:,p) = rotMat(:,:,2,p).*lambda;
%     Late(:,:,p) = rotMat(:,:,3,p).*lambda;
%     After(:,:,p) = rotMat(:,:,4,p).*lambda;
end

%%
blockIdx = 1;
vectors = mean(rotMat.(block_name{blockIdx}),3);

figure(1); clf; hold on
plot([0 vectors(1,1)],[0 vectors(2,1)],'LineWidth',1.5)
plot([0 vectors(1,2)],[0 vectors(2,2)],'LineWidth',1.5)
plot([0 1],[0 0],'k')
plot([0 0],[0 1],'k')
axis([-0.45 1 -0.45 1])
axis square

%%
subj = 4;
blockIdx = 4;

plot_subj(data{subj}.(block_name{blockIdx}),0)
subplot(2,2,1)
plot(rotMat(1,1,blockIdx,subj,group)*phasor(:,subj,group),'k','LineWidth',3)
subplot(2,2,2)
plot(rotMat(1,2,blockIdx,subj,group)*phasor(:,subj,group),'k','LineWidth',3)
subplot(2,2,3)
plot(rotMat(2,1,blockIdx,subj,group)*phasor(:,subj,group),'k','LineWidth',3)
subplot(2,2,4)
plot(rotMat(2,2,blockIdx,subj,group)*phasor(:,subj,group),'k','LineWidth',3)

%%
blockIdx = 1;
mat = mean(rotMat.(block_name{blockIdx}),3);

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
imagesc(mat,clims)
colormap(map)
set(gca,'TickDir','out','Xtick',[],'Ytick',[])
axis square

%%
function e = scale(params,cplx_data,Nfreq)
    cplx_fit = params(1:Nfreq);
    theta = [1; params(Nfreq+1:end)];
    Nblocks = length(theta)/4;
    theta = reshape(theta,[4 Nblocks]);
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