% Performs analysis and makes figures for the gain matrices. The argument
% "experiment" changes how figures are plotted based on whether "data" is
% from experiment 1 or 2. 
% 
% This function also generates "gain_matrix.csv" for experiment 1 and
% "gain_matrix2.csv" for experiment 2. These files contain the off-diagonal
% values of the gain matrices for statistical analysis in R. To generate
% these files, uncomment lines _____

function graph_gainMatrix(data,experiment)

% set variables for analysis
groups = {'rot','mir'};
blocks = {'baseline','pert1','pert2','pert3','pert4','post'};
graph_name = {'Baseline','Early','Train2','Train3','Late','Post'};
names = {'x_x','x_y','y_x','y_y'};
Nblock = length(blocks);
Nsubj = length(data.rot);
Nfreq = length(data.rot{1}.(blocks{1}).freqX);
Ngroup = length(groups);
Ntrials = size(data.rot{1}.(blocks{1}).Rhand.x_pos,2);

paramsInit = zeros([2*Nblock 1]); % initialize parameters as zeros

% this loop computes the gain matrices
for q = 1:Ngroup % loop over groups
    for p = 1:Nsubj % loop over subjects
        for k = 1:Nblock % loop over blocks
            % loop over combinations of target and hand axes
            % i=1: x-target to x-hand
            % i=2: x-target to y-hand
            % i=3: y-target to x-hand
            % i=4: y-target to y-hand
            for i = 1:4
                % store complex ratios in cplx
                cplx(:,:,i,k,p,q) = data.(groups{q}){p}.(blocks{k})...
                    .Rhand.phasors.(names{i}).ratio;
            end
        end
        
        for i = 1:Nfreq
            
            % compute average baseline phasor for a particular frequency
            dat = mean(cplx(i,:,[1 4],1,p,q),2);
            
            % x-component of xx and yy baseline phasor
            x = cos(angle(dat));
            
            % y-component of xx and yy baseline phasor
            y = sin(angle(dat));
            
            % project data onto baseline x phasor
            dat = cplx(i,:,[1 2],:,p,q); % extract data for x-target to x-/y-hand 
            phasor = reshape(dat, [1 numel(dat)]);
            num = [real(phasor); imag(phasor)]; % separate real and imaginary parts of phasors
            unitVec = [x(1); y(1)]; % unit vector with phase equal to baseline vector
            gainX(:,i,p,q) = dot(num, repmat(unitVec, [1 numel(dat)])); % project num onto unitVec
            frac1(:,i,p,q) = abs(gainX(:,i,p,q))./abs(phasor)'; % calculate proportion of gain retained in projection
            
            % project data onto baseline y phasor
            dat = cplx(i,:,[3 4],:,p,q); % extract data for x-target to x-/y-hand 
            phasor = reshape(dat, [1 numel(dat)]);
            num = [real(phasor); imag(phasor)]; % separate real and imaginary parts of phasors
            unitVec = [x(2); y(2)]; % unit vector with phase equal to baseline vector
            gainY(:,i,p,q) = dot(num, repmat(unitVec, [1 numel(dat)])); % project num onto unitVec
            frac2(:,i,p,q) = abs(gainY(:,i,p,q))./abs(phasor)'; % calculate proportion of gain retained in projection
        end
    end
end

% combine opt1 and opt2
thetaOpt = [reshape(gainX,[Ntrials 2 Nblock Nfreq Nsubj Ngroup]) ...
    reshape(gainY,[Ntrials 2 Nblock Nfreq Nsubj Ngroup])]; 
thetaOpt = permute(thetaOpt, [2 3 4 1 5 6]);

% if desired, weight estimated gains by the amount of gain lost in the
% projection
% fracOpt = [reshape(frac1,[Ntrials 2 Nblock Nfreq Nsubj Ngroup]) ...
%     reshape(frac2,[Ntrials 2 Nblock Nfreq Nsubj Ngroup])];
% fracOpt = permute(fracOpt, [2 3 4 1 5 6]);
% thetaOpt = thetaOpt .* (1./fracOpt);

% shape thetaOpt into gain matrix format: the first and second dimensions
% of gainMat correspond to the gain matrix for a given frequency (third 
% dimension), block (fourth dimension), subject (fifth dimension), and 
% group (sixth dimension)
gainMat = reshape(thetaOpt,[2 2 Nblock Nfreq Ntrials Nsubj Ngroup]); 
gainMat = permute(gainMat,[1 2 4 3 5 6 7]);

% fit theta to rotation matrices
for p = 1:Nsubj
    for m = 1:Ntrials
        for k = 1:Nblock
            for i = 1:Nfreq
                H = eye(2)*gainMat(:,:,i,k,m,p,1)';
                [U,S,V] = svd(H);
                if det(H') >= 0 % if determinant >= 0, rotation matrix can be computed
                    R = V*U';
                    thetaFit(i,k,m,p) = atan2(R(2,1),R(1,1))*180/pi;
                else % if rotation matrix can't be computed, set value to NaN
                    thetaFit(i,k,m,p) = NaN;
                end
            end
        end
    end
end

% compute orthogonal gain for mirror-reversal group
R = rotz(-45); % 45 degree clockwise rotation matrix
R = R(1:2,1:2);
perpAxis = R*[1 0]';

for p = 1:Nsubj
    for m = 1:Ntrials
        for k = 1:Nblock
            for i = 1:Nfreq
                gainOrth(i,k,m,p) = perpAxis'*gainMat(:,:,i,k,m,p,2)*perpAxis;
            end
        end
    end
end

thetaFit = reshape(permute(thetaFit,[3 2 1 4]),[Nblock*Ntrials Nfreq Nsubj]);
gainOrth = reshape(permute(gainOrth,[3 2 1 4]),[Nblock*Ntrials Nfreq Nsubj]);

%% plot vectors from gain matrices (Figure 5A, Figure 5-supplement 2A, Figure 6B, and Figure 6-supplement 1A)
% generate color maps
col1 = [0 128 0]/255;
col2 = [152 251 152]/255;
map1 = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2)...
    ,Nfreq)', linspace(col1(3),col2(3),Nfreq)'];

col1 = [128 0 128]/255;
col2 = [230 230 250]/255;
map2 = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2)...
    ,Nfreq)', linspace(col1(3),col2(3),Nfreq)'];

gblocks = [1 2 5 6]; % blocks to plot
trial = 1; % trial to plot

if experiment == 1 % plots Figure 5A and Figure 5-supplement 2A
    f = 14;
    subj = [7 9];
else % plots Figure 6B and Figure 6-supplement 1A
    f = 22;
    subj = [3 3];
end

% This loop plots average across subjects if p=1 and single subject data if
% p=2
for p = 1:2
    figure(f+p); clf
    if p == 1 % average data within group
        mat = squeeze(mean(thetaOpt(:,:,:,trial,:,:),5));
    end
    for i = 1:length(gblocks) % iterate over the blocks to plot
        if p == 2 % make plot for single subject in rotation group
            mat = squeeze(thetaOpt(:,:,:,trial,subj(1),:));
        end
        subplot(2,length(gblocks),i); hold on
        plot([0 1],[0 0],'k') % unit x vector
        plot([0 0],[0 1],'k') % unit y vector
        for k = 1:Nfreq
            plot([0 mat(1,gblocks(i),k,1)],[0 mat(2,gblocks(i),k,1)]...
                ,'LineWidth',1,'Color',map1(k,:)) % plot green vectors
            plot([0 mat(3,gblocks(i),k,1)],[0 mat(4,gblocks(i),k,1)]...
                ,'LineWidth',1,'Color',map2(k,:)) % plot purple vectors
        end
        axis([-0.95 1.5 -0.95 1.5])
        axis square
        title(graph_name(gblocks(i)))
        if i == 1
            ylabel('Rotation')
        end
        
        if p == 2 % make plot for single subject in mirror-reversal group
            mat = squeeze(thetaOpt(:,:,:,trial,subj(2),:));
        end
        subplot(2,length(gblocks),i+length(gblocks)); hold on
        plot([0 1],[0 0],'k') % unit x vector
        plot([0 0],[0 1],'k') % unit y vector
        for k = 1:Nfreq
            plot([0 mat(1,gblocks(i),k,2)],[0 mat(2,gblocks(i),k,2)]...
                ,'LineWidth',1,'Color',map1(k,:)) % plot green vectors
            plot([0 mat(3,gblocks(i),k,2)],[0 mat(4,gblocks(i),k,2)]...
                ,'LineWidth',1,'Color',map2(k,:)) % plot purple vectors
        end
        axis([-0.95 1.5 -0.95 1.5])
        axis square
        if i == 1
            ylabel('Mirror-Reversal')
        end
    end
end

%% plot off-diagonal elements (Figure 5B, Figure 5-supplement 2B, Figure 6C, and Figure 6-supplement 1B)
col = copper;
colors = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

% average the off-diagonal elements together
mat2 = permute(thetaOpt,[1 4 2 3 5 6]);
mat2 = reshape(mat2,[4 Nblock*Ntrials Nfreq Nsubj Ngroup]);
mat2(3,:,:,:,1) = -mat2(3,:,:,:,1); % for the rotation group, flip the sign of element 1,2
mat2 = [mean(mat2([1 4],:,:,:,:),1); mean(mat2(2:3,:,:,:,:),1)]; % average on-diagonal and off-diagonal elements
on = squeeze(nanmean(mat2(1,:,:,:,:),4)); % extract on-diagonal mean
off = squeeze(nanmean(mat2(2,:,:,:,:),4)); % extract off-diagonal mean

% This section generates csv files for statistical analysis in R. Uncomment
% the lines below to generate these matrices
% if experiment == 1
%     z = permute(mat2(2,[1 40 41],:,:,:),[4 3 2 5 1]);
%     z = reshape(z,[numel(z) 1]);
%     dlmwrite('gain_matrix.csv',z);
% elseif experiment == 2
%     z = permute(mat2(2,[1 30 31],:,:,:),[4 3 2 5 1]);
%     z = reshape(z,[numel(z) 1]);
%     dlmwrite('gain_matrix2.csv',z);
% end

% compute standard error over the averaged on- and off-diagonal elements
onAll = permute(mat2(1,:,:,:,:),[2 4 3 5 1]);
offAll = permute(mat2(2,:,:,:,:),[2 4 3 5 1]);
onSE = squeeze(nanstd(onAll,[],2)./sqrt(Nsubj-sum(isnan(onAll),2)));
offSE = squeeze(nanstd(offAll,[],2)./sqrt(Nsubj-sum(isnan(offAll),2)));

gblocks = Nblock*Ntrials; % blocks to plot
if experiment == 1 % plot Figure 5B and Figure 5-supplement 2B
    f = 17;
    ax = [1 gblocks -0.15 0.75];
    subj = [7 9];
else % plot Figure 6C and Figure 6-supplement 1B
    f = 25;
    ax = [1 gblocks -0.25 0.75];
    subj = [3 3];
end

% plot data averaged across subjects
figure(f); clf
for i = 1:Nfreq
    % plot rotation data
    subplot(1,2,1); hold on
    if i == 1 % plot background line and rectangle behind data
        plot([1 gblocks],[0 0],'k')
        rectangle('Position',[Ntrials+1 -0.5 4*Ntrials-1 2],'FaceColor',...
            [0 0 0 0.1],'EdgeColor','none');
    end
    for k = 1:Nblock % plot data from each block and frequency as different lines
        plotIdx = Ntrials*(k-1)+1:Ntrials*(k-1)+Ntrials;
        s = shadedErrorBar(plotIdx,off(plotIdx,i,1),offSE(plotIdx,i,1));
        editErrorBar(s,colors(i,:),1.5);
    end
    set(gca,'TickDir','out')
    axis(ax)
    yticks(0:0.3:0.6)
    xticks(1:Ntrials:41)
    xlabel('Trial Number')
    title('Rotation')
    if i == 1
        ylabel('Off-diagonal gain')
    end
    
    % plot mirror-reversal data
    subplot(1,2,2); hold on
    if i == 1 % plot background line and rectangle behind data
        plot([1 gblocks],[0 0],'k')
        rectangle('Position',[Ntrials+1 -0.5 4*Ntrials-1 2],'FaceColor',...
            [0 0 0 0.1],'EdgeColor','none');
    end
    for k = 1:Nblock % plot data from each block and frequency as different lines
        plotIdx = Ntrials*(k-1)+1:Ntrials*(k-1)+Ntrials;
        s = shadedErrorBar(plotIdx,off(plotIdx,i,2),offSE(plotIdx,i,2));
        editErrorBar(s,colors(i,:),1.5);
    end
    set(gca,'TickDir','out')
    axis(ax)
    yticks(0:0.3:0.6)
    xticks(1:Ntrials:41)
    xlabel('Trial Number')
    title('Mirror Reversal')
end

% plot single subject data
off(:,:,1) = squeeze(mat2(2,:,:,subj(1),1));
off(:,:,2) = squeeze(mat2(2,:,:,subj(2),2));

figure(f+1); clf
for i = 1:Nfreq
    % plot rotation data
    subplot(1,2,1); hold on
    if i == 1 % plot background line and rectangle behind data
        plot([1 gblocks],[0 0],'k')
        rectangle('Position',[Ntrials+1 -0.5 4*Ntrials-1 2],'FaceColor',...
            [0 0 0 0.1],'EdgeColor','none');
    end
    for k = 1:Nblock % plot data from each block and frequency as different lines
        plotIdx = Ntrials*(k-1)+1:Ntrials*(k-1)+Ntrials;
        plot(plotIdx,off(plotIdx,i,1),'Color',colors(i,:),'LineWidth',1.5);
    end
    set(gca,'TickDir','out')
    axis([1 gblocks -0.25 0.85])
    yticks(0:0.3:0.6)
    xticks(1:Ntrials:41)
    xlabel('Trial Number')
    title('Rotation')
    if i == 1
        ylabel('Off-diagonal gain')
    end
    
    % plot mirror-reversal data
    subplot(1,2,2); hold on
    if i == 1 % plot background line and rectangle behind data
        plot([1 gblocks],[0 0],'k')
        rectangle('Position',[Ntrials+1 -0.5 4*Ntrials-1 2],'FaceColor',...
            [0 0 0 0.1],'EdgeColor','none');
    end
    for k = 1:Nblock % plot data from each block and frequency as different lines
        plotIdx = Ntrials*(k-1)+1:Ntrials*(k-1)+Ntrials;
        plot(plotIdx,off(plotIdx,i,2),'Color',colors(i,:),'LineWidth',1.5);
    end
    set(gca,'TickDir','out')
    axis([1 gblocks -0.25 0.85])
    yticks(0:0.3:0.6)
    xticks(1:Ntrials:41)
    xlabel('Trial Number')
    title('Mirror Reversal')
end

%% plot rotation angle and orthogonal gain (Figure 5C, Figure 5-supplement 2C, Figure 6D, and Figure 6-supplement 1C)
col = copper;
colors = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

% compute mean and standard deviation of the rotation angle/orthogonal gain
gainOrthMu = mean(gainOrth,3);
gainOrthSE = std(gainOrth,[],3)/sqrt(Nsubj);
thetaFitMu = nanmean(thetaFit,3);
thetaFitSE = nanstd(thetaFit,[],3)./sqrt(Nsubj-sum(isnan(thetaFit),3));

if experiment == 1 % for Figure 5C and Figure 5-supplement 2C
    f = 19;
    ax1 = [1 Nblock*Ntrials -20 90];
    ax2 = [1 Nblock*Ntrials -1 1];
    ax3 = [1 Nblock*Ntrials -30 110];
    ax4 = [1 Nblock*Ntrials -1 1];
    subj = [7 9];
else % for Figure 6D and Figure 6-supplement 1C
    f = 27;
    ax1 = [1 Nblock*Ntrials -20 120];
    ax2 = [1 Nblock*Ntrials -1 1.15];
    ax3 = [1 Nblock*Ntrials -45 135];
    ax4 = [1 Nblock*Ntrials -1 1];
    subj = [3 3];
end

% plot data averaged across subjects
figure(f); clf

% rotation angle for rotation group
subplot(1,2,1); hold on
rectangle('Position',[Ntrials+1 -40 4*Ntrials-1 180],'FaceColor',[0 0 0 0.1],'EdgeColor','none');
plot([1 Nblock*Ntrials],[0 0],'--k','LineWidth',1) % ideal baseline response
plot([1 Nblock*Ntrials],[90 90],'--k','LineWidth',1) % ideal compensation
for k = 1:Nblock
    for i = 1:Nfreq
        plotIdx = Ntrials*(k-1)+1:Ntrials*(k-1)+Ntrials;
        s = shadedErrorBar(plotIdx,thetaFitMu(plotIdx,i),thetaFitSE(plotIdx,i));
        editErrorBar(s,colors(i,:),1.5);
    end
end
set(gca,'TickDir','out')
title('Rotation')
xticks(1:Ntrials:41)
yticks(-30:30:120)
ylabel(['Angle (' char(176) ')'])
axis(ax1)

% orthogonal gain for mirror-reversal group
subplot(1,2,2); hold on
rectangle('Position',[Ntrials+1 -2 4*Ntrials-1 4],'FaceColor',[0 0 0 0.1],'EdgeColor','none');
plot([1 Nblock*Ntrials],[1 1],'--k','LineWidth',1,'HandleVisibility','off') % ideal baseline response
plot([1 Nblock*Ntrials],[-1 -1],'--k','LineWidth',1,'HandleVisibility','off') % ideal compensation
plot([1 Nblock*Ntrials],[0 0],'k','HandleVisibility','off')
for k = 1:Nblock
    for i = 1:Nfreq
        plotIdx = Ntrials*(k-1)+1:Ntrials*(k-1)+Ntrials;
        s = shadedErrorBar(plotIdx,gainOrthMu(plotIdx,i),gainOrthSE(plotIdx,i));
        editErrorBar(s,colors(i,:),1.5);
    end
end
set(gca,'TickDir','out')
title('Mirror-Reversal')
xticks(1:Ntrials:41)
yticks(-1:0.5:1)
ylabel('Gain orthogonal to mirroring axis')
axis(ax2)

% plot single subject data
thetaFitMu = thetaFit(:,:,subj(1));
gainOrthMu = gainOrth(:,:,subj(2));

figure(f+1); clf

% rotation angle for rotation group
subplot(1,2,1); hold on
rectangle('Position',[Ntrials+1 -180 4*Ntrials-1 360],'FaceColor',[0 0 0 0.1],'EdgeColor','none');
plot([1 Nblock*Ntrials],[0 0],'--k','LineWidth',1) % ideal baseline response
plot([1 Nblock*Ntrials],[90 90],'--k','LineWidth',1) % ideal compensation
for k = 1:Nblock
    for i = 1:Nfreq
        plotIdx = Ntrials*(k-1)+1:Ntrials*(k-1)+Ntrials;
        plot(plotIdx,thetaFitMu(plotIdx,i),'Color',colors(i,:),'LineWidth',1.5);
    end
end
set(gca,'TickDir','out')
title('Rotation')
xticks(1:Ntrials:41)
yticks(-180:45:180)
ylabel(['Angle (' char(176) ')'])
axis(ax3)

% orthogonal gain for mirror-reversal group
subplot(1,2,2); hold on
rectangle('Position',[Ntrials+1 -2 4*Ntrials-1 4],'FaceColor',[0 0 0 0.1],'EdgeColor','none');
plot([1 Nblock*Ntrials],[1 1],'--k','LineWidth',1,'HandleVisibility','off') % ideal baseline response
plot([1 Nblock*Ntrials],[-1 -1],'--k','LineWidth',1,'HandleVisibility','off') % ideal compensation
plot([1 Nblock*Ntrials],[0 0],'k','HandleVisibility','off')
for k = 1:Nblock
    for i = 1:Nfreq
        plotIdx = Ntrials*(k-1)+1:Ntrials*(k-1)+Ntrials;
        plot(plotIdx,gainOrthMu(plotIdx,i),'Color',colors(i,:),'LineWidth',1.5);
    end
end
set(gca,'TickDir','out')
title('Mirror-Reversal')
xticks(1:Ntrials:41)
yticks(-1:0.5:1)
ylabel('Gain orthogonal to mirroring axis')
axis(ax4)

%% plot gain matrices (Figure 5-supplement 1)
if experiment == 1 % only plot the gain matrices for the first experiment

    % generate color maps for the matrices
    col1 = [1 0 0];
    col2 = [1 1 1];
    Nstep = 100;
    map1 = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2)...
        ,Nstep)', linspace(col1(3),col2(3),Nstep)'];
    
    col1 = [1 1 1];
    col2 = [0 0 1];
    map2 = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2)...
        ,Nstep)', linspace(col1(3),col2(3),Nstep)'];
    
    map = [map1; map2];
    clims = [-1 1];
    
    % choose whether to make plot for single subject...
    % subj = 6;
    % mat = squeeze(thetaOpt(:,:,:,subj,:));
    
    % OR average across subjects
    mat = squeeze(mean(gainMat,6));
    
    % choose which trial to plot
    trial = 1;
    
    % i=1 plots rotation data, i=2 plots mirror-reversal data
    for i = 1:Ngroup
        figure(i+20); clf
        for k = 1:Nblock
            for p = 1:Nfreq
                subplot(Nblock,Nfreq,(k-1)*Nfreq+p)
                imagesc(mat(:,:,p,k,trial,i),clims)
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
end