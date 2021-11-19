% Performs analysis and makes figures for the gain matrices. The argument
% "experiment" changes how figures are plotted based on whether "data" is
% from experiment 1 or 2. 
% 
% This function also generates "gain_matrix.csv" for experiment 1 and
% "gain_matrix2.csv" for experiment 2. These files contain the off-diagonal
% values of the gain matrices for statistical analysis in R. To generate
% these files, uncomment lines _____

function graph_gainMatrix(data)

% set variables for analysis
groups = {'rot','mir'};
blocks = {'B1','B2','B3'};
graph_name = {'Baseline','Early','Late'};
names = {'Hur','Hud','B','F'};
Nblock = length(blocks);
% Nsubj = length(data.rot);
Nfreq = length(data.rot{1}.(blocks{1}).freq);
Ngroup = length(groups);
Ntrial = 4;

for i = 1:Ngroup
    for j = 1:length(names)
        cplx.(groups{i}).(names{j}) = NaN(4,Nfreq,Ntrial,Nblock,length(data.(groups{i})));
    end
end
% cplx.Hud = NaN(4,Nfreq,Ntrial,Nblock,Nsubj,Ngroup);
% cplx.B = NaN(4,Nfreq,Ntrial,Nblock,Nsubj,Ngroup);
% cplx.F = NaN(4,Nfreq,Ntrial,Nblock,Nsubj,Ngroup);

% this loop computes the gain matrices
for q = 1:Ngroup % loop over groups
    Nsubj = length(data.(groups{q}));
    for p = 1:Nsubj % loop over subjects
        for k = 1:Nblock % loop over blocks
            % loop over combinations of target and hand axes
            % i=1: x-target to x-hand
            % i=2: x-target to y-hand
            % i=3: y-target to x-hand
            % i=4: y-target to y-hand
            
            a = data.(groups{q}){p}.(blocks{k});
            n = size(a.B,4);
            
            for j = 1:4
                cplx.(groups{q}).(names{j})(:,:,1:n,k,p) = reshape(a.(names{j}),[4 Nfreq n]);
            end
%             Hur(:,:,1:n,k,p,q) = reshape(a.Hur,[4 Nfreq n]);
%             Hud(:,:,1:n,k,p,q) = reshape(a.Hud,[4 Nfreq n]);
%             F(:,:,1:n,k,p,q) = reshape(a.F,[4 Nfreq n]);
%             B(:,:,1:n,k,p,q) = reshape(a.B,[4 Nfreq n]);
        end
        
        for j = 1:4
            for i = 1:Nfreq

                % compute average baseline phasor for a particular frequency
                dat = mean(cplx.(groups{q}).(names{j})([1 4],i,:,1,p),3,'omitnan');

                % x-component of xx and yy baseline phasor
                x = cos(angle(dat));

                % y-component of xx and yy baseline phasor
                y = sin(angle(dat));

                % project data onto baseline x phasor
                phasor = cplx.(groups{q}).(names{j})([1 2],i,:,:,p); % extract data for x-target to x-/y-hand 
                phasor = reshape(phasor,[1 numel(phasor)]);
                num = [real(phasor); imag(phasor)]; % separate real and imaginary parts of phasors
                unitVec = [x(1); y(1)]; % unit vector with phase equal to baseline vector
                gainX.(groups{q}).(names{j})(:,i,p) = dot(num, repmat(unitVec, [1 numel(phasor)])); % project num onto unitVec
    %             frac1(:,i,p,q) = abs(gainX(:,i,p,q))./abs(phasor)'; % calculate proportion of gain retained in projection

                % project data onto baseline y phasor
                phasor = cplx.(groups{q}).(names{j})([3 4],i,:,:,p); % extract data for x-target to x-/y-hand 
                phasor = reshape(phasor,[1 numel(phasor)]);
                num = [real(phasor); imag(phasor)]; % separate real and imaginary parts of phasors
                unitVec = [x(2); y(2)]; % unit vector with phase equal to baseline vector
                gainY.(groups{q}).(names{j})(:,i,p) = dot(num, repmat(unitVec, [1 numel(phasor)])); % project num onto unitVec
    %             frac2(:,i,p,q) = abs(gainY(:,i,p,q))./abs(phasor)'; % calculate proportion of gain retained in projection
            end
        end
    end
end

% compute orthogonal gain for mirror-reversal group
R = rotz(-45); % 45 degree clockwise rotation matrix
R = R(1:2,1:2);
perpAxis = R*[1 0]';

% combine opt1 and opt2
for z = 1:Ngroup
    Nsubj = length(data.(groups{z}));
    for j = 1:4
        thetaOpt.(groups{z}).(names{j}) = [reshape(gainX.(groups{z}).(names{j}),[2 4 Nblock Nfreq Nsubj]); ...
            reshape(gainY.(groups{z}).(names{j}),[2 4 Nblock Nfreq Nsubj])];
        thetaOpt.(groups{z}).(names{j}) = permute(thetaOpt.(groups{z}).(names{j}), [1 3 4 2 5]);
        
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
        gainMat.(groups{z}).(names{j}) = reshape(thetaOpt.(groups{z}).(names{j}),[2 2 Nblock Nfreq Ntrial Nsubj]);
        gainMat.(groups{z}).(names{j}) = permute(gainMat.(groups{z}).(names{j}),[1 2 4 3 5 6]);
    end
end

for j = 1:4        
    % fit theta to rotation matrices
    for p = 1:Nsubj
        for m = 1:Ntrial
            for k = 1:Nblock
                for i = 1:Nfreq
                    H = gainMat.rot.(names{j})(:,:,i,k,m,p)';
                    if sum(isnan(H),'all') == 0
                        [U,S,V] = svd(H);
                        if det(H') >= 0 % if determinant >= 0, rotation matrix can be computed
                            R = V*U';
                            thetaFit.(names{j})(i,k,m,p) = atan2(R(2,1),R(1,1))*180/pi;
                        else % if rotation matrix can't be computed, set value to NaN
                            thetaFit.(names{j})(i,k,m,p) = NaN;
                        end
                    else
                        thetaFit.(names{j})(i,k,m,p) = NaN;
                    end
                end
            end
        end
    end
    
    for p = 1:Nsubj
        for m = 1:Ntrial
            for k = 1:Nblock
                for i = 1:Nfreq
                    gainOrth.(names{j})(i,k,m,p) = perpAxis'*gainMat.(groups{z}).(names{j})(:,:,i,k,m,p)*perpAxis;
                end
            end
        end
    end
    
    thetaFit.(names{j}) = reshape(permute(thetaFit.(names{j}),[3 2 1 4]),[Nblock*Ntrial Nfreq Nsubj]);
    gainOrth.(names{j}) = reshape(permute(gainOrth.(names{j}),[3 2 1 4]),[Nblock*Ntrial Nfreq Nsubj]);
end

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

trial = 1; % trial to plot

for j = 1:4
    figure(j); clf
    for q = 1:Ngroup
        mat = squeeze(mean(thetaOpt.(groups{q}).(names{j})(:,:,:,trial,:),5));
        
        for i = 1:Nblock % iterate over the blocks to plot
            subplot(2,Nblock,(q-1)*Nblock+i); hold on
            plot([0 1],[0 0],'k') % unit x vector
            plot([0 0],[0 1],'k') % unit y vector
            for k = 1:2:Nfreq
                plot([0 mat(1,i,k)],[0 mat(2,i,k)]...
                    ,'LineWidth',1,'Color',map1(k,:)) % plot green vectors
                plot([0 mat(3,i,k)],[0 mat(4,i,k)]...
                    ,'LineWidth',1,'Color',map2(k,:)) % plot purple vectors
            end
            axis([-0.95 1.5 -0.95 1.5])
            axis square
            title(graph_name(i))
            if i == 1
                if q == 1
                    ylabel('Rotation')
                else
                    ylabel('Mirror')
                end
            end
            
%             subplot(2,Nblock,Nblock+i); hold on
%             plot([0 1],[0 0],'k') % unit x vector
%             plot([0 0],[0 1],'k') % unit y vector
%             for k = 1:Nfreq
%                 plot([0 mat(1,i,k)],[0 mat(2,i,k)]...
%                     ,'LineWidth',1,'Color',map1(k,:)) % plot green vectors
%                 plot([0 mat(3,i,k)],[0 mat(4,i,k)]...
%                     ,'LineWidth',1,'Color',map2(k,:)) % plot purple vectors
%             end
%             axis([-0.95 1.5 -0.95 1.5])
%             axis square
%             if i == 1
%                 ylabel('Mirror-Reversal')
%             end
        end
    end
end

%% plot rotation angle and orthogonal gain (Figure 5C, Figure 5-supplement 2C, Figure 6D, and Figure 6-supplement 1C)
col = copper;
colors = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

figure(5); clf
for j = 1:4
    % compute mean and standard deviation of the rotation angle/orthogonal gain
    gainOrthMu = mean(gainOrth.(names{j}),3);
    gainOrthSE = std(gainOrth.(names{j}),[],3)/sqrt(Nsubj);
    thetaFitMu = nanmean(thetaFit.(names{j}),3);
    thetaFitSE = nanstd(thetaFit.(names{j}),[],3)./sqrt(Nsubj-sum(isnan(thetaFit.(names{j})),3));
    
    % rotation angle for rotation group
    subplot(2,4,j); hold on
    plot([1 Nblock*Ntrial],[0 0],'--k','LineWidth',1) % ideal baseline response
    plot([1 Nblock*Ntrial],[90 90],'--k','LineWidth',1) % ideal compensation
    for k = 1:Nblock
        for i = 1:Nfreq
            plotIdx = Ntrial*(k-1)+1:Ntrial*(k-1)+Ntrial;
            s = shadedErrorBar(plotIdx,thetaFitMu(plotIdx,i),...
                thetaFitSE(plotIdx,i));
            editErrorBar(s,colors(i,:),1.5);
        end
    end
    set(gca,'TickDir','out')
    title(names{j})
    xticks(1:Ntrial:9)
    ylim([-50 120])
%     yticks(-30:30:120)
    if j == 1
        ylabel(['Angle (' char(176) ')'])
    end

    % orthogonal gain for mirror-reversal group
    subplot(2,4,j+4); hold on
    plot([1 Nblock*Ntrial],[1 1],'--k','LineWidth',1,'HandleVisibility','off') % ideal baseline response
    plot([1 Nblock*Ntrial],[-1 -1],'--k','LineWidth',1,'HandleVisibility',...
        'off') % ideal compensation
    plot([1 Nblock*Ntrial],[0 0],'k','HandleVisibility','off')
    for k = 1:Nblock
        for i = 1:Nfreq
            plotIdx = Ntrial*(k-1)+1:Ntrial*(k-1)+Ntrial;
            s = shadedErrorBar(plotIdx,gainOrthMu(plotIdx,i),...
                gainOrthSE(plotIdx,i));
            editErrorBar(s,colors(i,:),1.5);
        end
    end
    set(gca,'TickDir','out')
    xticks(1:Ntrial:9)
    ylim([-1.2 1.2])
%     yticks(-1:0.5:1)
    if j == 1
        ylabel('Gain orthogonal to mirroring axis')
    end
end

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