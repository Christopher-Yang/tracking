output = 'Lhand';
groups = {'rot','mir'};
names = {'x_x_all','x_y_all','y_x_all','y_y_all'};
Nblock = length(block_name);
Nsubj = length(data)-1;
Nfreq = length(data{end}.freqX);
Ngroups = 1;
paramsInit = zeros([2 Nblock]);
idx = find(contains(graph_name,'('));

for p = 1:Nsubj % loop over subjects
    for k = 1:Nblock
        for i = 1:4
            cplx(:,i,k,p) = mean(data{p}.(block_name{k}).phasors.(output).(names{i}).ratio,2); % put all complex ratios in cplx
        end
    end
    
    for i = 1:Nfreq
        x = cos(angle(cplx(i,[1 4],1,p)));
        y = sin(angle(cplx(i,[1 4],1,p)));
        template(:,i,p) = x+y*1j; % generate template of gain = 1 and phase same as baseline
        error = @(params) scale(params,cplx(i,[1 2],:,p),template(1,i,p));
        opt1(:,:,i,p) = fmincon(error,paramsInit); % fit gains for x frequencies
        error = @(params) scale(params,cplx(i,[3 4],:,p),template(2,i,p));
        opt2(:,:,i,p) = fmincon(error,paramsInit); % fit gains for y frequencies
    end
end

% combine opt1 and opt2
thetaOpt = [reshape(opt1,[2 Nblock Nfreq Nsubj Ngroups]); reshape(opt2,[2 Nblock Nfreq Nsubj Ngroups])];

% shape thetaOpt into gain matrix format
rotMat = reshape(thetaOpt,[2 2 Nblock Nfreq Nsubj Ngroups]);
rotMat = permute(rotMat,[1 2 4 3 5 6]);

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

%% display results of fitting process
g = 1;
subj = 1;
blockIdx = 1;
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
plot_subj(data{subj}.(block_name{blockIdx}).phasors.(output))
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
% for averaging across subjects
m = mean(thetaOpt,4);
for i = 1:Nfreq
    mat(:,:,i) = squeeze(m(:,:,i));
end

if strcmp(output,'cursor')
    names2 = {'X_{T}X_{O} (response)','X_{T}Y_{O}','Y_{T}X_{O}','Y_{T}Y_{O} (response)'};
elseif strcmp(output,'Rhand')
    names2 = {'X_{T}X_{O}','X_{T}Y_{O}','Y_{T}X_{O} (response)','Y_{T}Y_{O}'};
elseif strcmp(output,'Lhand')
    names2 = {'X_{T}X_{O}','X_{T}Y_{O} (response)','Y_{T}X_{O}','Y_{T}Y_{O}'};
end

index = 1:Nblock;
remove = 1;
if remove
    index(idx) = [];
end

figure(2); clf
for k = 1:Nfreq
    for i = index
        subplot(1,Nfreq,k); hold on
        plot([0 mat(1,i,k)],[0 mat(3,i,k)],'LineWidth',1.5,'Color',map1(i,:))
        plot([0 mat(2,i,k)],[0 mat(4,i,k)],'LineWidth',1.5,'Color',map2(i,:))
        plot([0 1],[0 0],'k')
        plot([0 0],[0 1],'k')
        axis([-0.45 1 -0.45 1])
        axis square
    end
    if k == 1
        title('Low Freq')
    elseif k == 6
        title('High Freq')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gain matrices
figure(3); clf
for i = 1:Nfreq
    subplot(1,Nfreq,i)
    imagesc(mat(:,index,i),clims)
    colormap(map)
    set(gca,'TickDir','out','Xtick',[],'Ytick',[])
    if i == 1
        title('Low freq')
        yticks(1:4)
        yticklabels(names2)
    elseif i == Nfreq
        title('High freq')
        if remove == 0
            xticks(idx)
            xticklabels(graph_name(idx))
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
    phasors = params*template;
    e = (phasors - squeeze(cplx_data)).^2;
    for k = 1:numel(e)
        e(k) = norm(e(k));
    end
    e = sum(sum(e));
end