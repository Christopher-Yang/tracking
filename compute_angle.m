Nsubj = length(data)-1;
Nblock = length(block_name);
names = {'x_x_all','x_y_all','y_x_all','y_y_all'};
output = 'Rhand';

for p = 1:Nsubj
    for k = 1:length(block_name)
        for i = 1:length(names)
            cplx(:,i,k) = mean(data{p}.(block_name{k}).phasors.(output).(names{i}).ratio,2);
        end
    end
    
    theta = zeros([4*length(block_name)-1 1]);
    paramsInit = [cplx(:,1,1); theta];
    
    error = @(params) scale(params,cplx);
    paramsOpt = fmincon(error,paramsInit);
    phasor(:,p) = paramsOpt(1:6);
    thetaOpt = [1; paramsOpt(7:end)];
    thetaOpt = reshape(thetaOpt, [4 length(block_name)]);
    
    for k = 1:length(block_name)
        rotMat(:,:,k,p) = reshape(thetaOpt(:,k), [2 2]);
    end
    
    ph = abs(phasor(:,p));
    lambda = mean(ph);
    
    gainMat(:,:,:,p) = rotMat(:,:,:,p).*lambda;
end

% for plotting heatmaps
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
col1 = [255 192 203]/255;
col2 = [139 0 0]/255;
map1 = [linspace(col1(1),col2(1),Nblock)', linspace(col1(2),col2(2),Nblock)', linspace(col1(3),col2(3),Nblock)'];

col1 = [176 224 230]/255;
col2 = [0 0 205]/255;
map2 = [linspace(col1(1),col2(1),Nblock)', linspace(col1(2),col2(2),Nblock)', linspace(col1(3),col2(3),Nblock)'];

%% plot gain matrices and vectors
% subj = 2;
% block = 1:4;

% for single subjects
% mat = gainMat(:,:,block,subj);

% for averaging across subjects
% mat = mean(gainMat(:,:,block,:),4);

m = mean(gainMat,4);
mat(1,:) = m(1,1,:);
mat(2,:) = m(1,2,:);
mat(3,:) = m(2,1,:);
mat(4,:) = m(2,2,:);

figure(1); clf
for i = 1:4
    % gain matrices
    subplot(4,1,i)
    imagesc(mat(i,:),clims)
    colormap(map)
    set(gca,'Xtick',[],'Ytick',[])
end

figure(2); clf; hold on
for i = 1:Nblock
    plot([0 m(1,1,i)],[0 m(2,1,i)],'LineWidth',1.5,'Color',map1(i,:))
    plot([0 m(1,2,i)],[0 m(2,2,i)],'LineWidth',1.5,'Color',map2(i,:))
    plot([0 1],[0 0],'k')
    plot([0 0],[0 1],'k')
    axis([-0.45 1 -0.45 1])
    axis square
end

%% fitted phasors for individual subjects
subj = 2;
blockIdx = 3;

plot_subj(data{subj}.(block_name{blockIdx}).phasors.(output),2)
subplot(2,2,1)
plot(rotMat(1,1,blockIdx,subj)*phasor(:,subj),'k','LineWidth',3)
subplot(2,2,2)
plot(rotMat(1,2,blockIdx,subj)*phasor(:,subj),'k','LineWidth',3)
subplot(2,2,3)
plot(rotMat(2,1,blockIdx,subj)*phasor(:,subj),'k','LineWidth',3)
subplot(2,2,4)
plot(rotMat(2,2,blockIdx,subj)*phasor(:,subj),'k','LineWidth',3)

%%
function e = scale(params,cplx_data)
    cplx_fit = params(1:6);
    theta = [1; params(7:end)];
    theta = reshape(theta,[4 length(theta)/4]);
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