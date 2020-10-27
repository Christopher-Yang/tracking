% function graph_gainMatrix(data)

graph_name = {'Baseline','Block 1','Block 2','Block 3','Block 4','Block 5'};

Nsubj = length(data);
Nblock = length(data{1});
Nfreq = length(data{1}{1}.sineParams.tX_freq{1});
Ntypes = 4;
paramsInit = zeros([2 Nblock]);

analysis = 2;

if analysis == 1
    names{1} = {'xTarg_x','xTarg_y','yTarg_x','yTarg_y'};
    names{2} = {'xCurs_x','xCurs_y','yCurs_x','yCurs_y'};
    names{3} = names{1};
    names{4} = names{2};
else
    names{1} = 'F';
    names{2} = 'B';
    names{3} = 'F2';
    names{4} = 'B2';
end

for q = 1:Ntypes
    for p = 1:Nsubj % loop over subjects
        for k = 1:Nblock % loop over block
            if k == 1
                idx = 1:4;
            else
                idx = [2 1 4 3];
            end
            
            if analysis == 1
                a = data{p}{k}.phasors;
                for i = 1:4 % loop over all combinations of hand/target axes
                    switch q
                        case 1
                            phasors = a.(names{q}{idx(i)}){1}.ratio;
                        case 2
                            phasors = a.(names{q}{idx(i)}){2}.ratio;
                        case 3
                            phasors = reshape(permute([a.(names{q}{idx(i)}){3}.ratio a.(names{q}{idx(i)}){4}.ratio],[2 1]),[Nfreq 1]);
                        case 4
                            phasors = reshape(permute([a.(names{q}{idx(i)}){4}.ratio a.(names{q}{idx(i)}){3}.ratio],[2 1]),[Nfreq 1]);
                    end
                    cplx(:,i,k,p,q) = phasors; % put all complex ratios in cplx
                end
            else
                a = reshape(data{p}{k}.(names{q}),[4 Nfreq]);
                for i = 1:4
                    cplx(:,i,k,p,q) = a(i,:);
                end
            end
        end
        
        for i = 1:Nfreq
            x = cos(angle(cplx(i,[1 4],1,p,q)));
            y = sin(angle(cplx(i,[1 4],1,p,q)));
            template(:,i,p,q) = x+y*1j; % generate template of gain = 1 and phase same as baseline
            error = @(params) scale(params,cplx(i,[1 2],:,p,q),template(1,i,p,q));
            opt1(:,:,i,p,q) = fmincon(error,paramsInit); % fit gains for x frequencies
            error = @(params) scale(params,cplx(i,[3 4],:,p,q),template(2,i,p,q));
            opt2(:,:,i,p,q) = fmincon(error,paramsInit); % fit gains for y frequencies
        end
    end
end

% combine opt1 and opt2
thetaOpt = [reshape(opt1,[2 Nblock Nfreq Nsubj Ntypes]); reshape(opt2,[2 Nblock Nfreq Nsubj Ntypes])];

% shape thetaOpt into gain matrix format
rotMat = reshape(thetaOpt,[2 2 Nblock Nfreq Nsubj Ntypes]);
rotMat = permute(rotMat,[1 2 4 3 5 6]);

% for averaging across subjects
mat = squeeze(mean(thetaOpt,4));

% for plotting lines
col = lines;
col = col(1:7,:);

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
col1 = [0 128 0]/255;
col2 = [152 251 152]/255;
map1 = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2),Nfreq)', linspace(col1(3),col2(3),Nfreq)'];

col1 = [128 0 128]/255;
col2 = [230 230 250]/255;
map2 = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2),Nfreq)', linspace(col1(3),col2(3),Nfreq)'];

%% plot vectors
if analysis == 1
    labels = {'Target sines (individual)','Cursor sines (individual)','Target sines (dual)','Cursor sines (dual)'};
else
    labels = {'Feedforward (indivudal)','Feedback (individual)','Feedforward (dual)','Feedback (dual)'};
end

figure(1); clf
for q = 1:Ntypes
    for k = 1:Nblock
        subplot(Ntypes,Nblock,Nblock*(q-1)+k); hold on
        plot([0 1],[0 0],'k')
        plot([0 0],[0 1],'k')
        for i = 1:Nfreq
            plot([0 mat(1,k,i,q)],[0 mat(2,k,i,q)],'LineWidth',1.5,'Color',map1(i,:))
            plot([0 mat(3,k,i,q)],[0 mat(4,k,i,q)],'LineWidth',1.5,'Color',map2(i,:))
            axis([-0.7 2 -0.7 2])
            axis square
        end
        xticks([])
        yticks([])
        
        if k == 1
            ylabel(labels{q})
            if q == 1
                title('Baseline')
            end
        else
            if q == 1
                title(['Block ' num2str(k-1)])
            end
        end
    end
end

%% plot gain matrices as 2x2
subj = 1;
type = 1;
rMat = rotMat(:,:,:,:,subj,type);

% rMat = mean(rotMat,5);

figure(4); clf
for i = 1:Nfreq
    for k = 1:Nblock
        subplot(Nfreq,Nblock,(i-1)*Nblock+k)
        imagesc(rMat(:,:,i,k),clims)
        colormap(map)
        pbaspect([1 1 1])
        set(gca,'Xtick',[],'Ytick',[])
        if i == 1
            title(graph_name(k))
        end
        if k == 1 && i == 1
            ylabel('Low freq')
        elseif k == 1 && i == Nfreq
            ylabel('High freq')
        end
    end
end

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

%% graph average fitted phasors across subjects
col1 = [1 0.9 0.3];
col2 = [1 0 0];
col = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2),Nfreq)', linspace(col1(3),col2(3),Nfreq)'];
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
        plot(dat_all(:,i,k),'.','Color',col(i,:),'MarkerSize',10)
        if i > 1
            plot(dat(i-1:i,k),'-o','Color',col(i,:),'LineWidth',1,'MarkerFaceColor',col(i-1,:),'MarkerEdgeColor','none')
        end
    end
    plot(phasorLines(:,k),'-ok','LineWidth',2,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none')
    plot([-1.5 1.5],[0 0],'k','LineWidth',1)
    plot([0 0],[-1.5 1.5],'k','LineWidth',1)
    axis([-1.25 1.25 -1.25 1.25])
    pbaspect([1 1 1])
end
% end

%%
function e = scale(params,cplx_data,template)
    phasors = params*template;
    e = (phasors - squeeze(cplx_data)).^2;
    for k = 1:numel(e)
        e(k) = norm(e(k));
    end
    e = sum(sum(e));
end