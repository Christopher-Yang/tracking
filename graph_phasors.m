Nsubj = length(data);
Nfreq = 12;
col = copper;
col = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

block = 'B3';
names1 = {'xTarg_x','yTarg_x','xTarg_y','yTarg_y'};
names2 = {'xCurs_x','yCurs_x','xCurs_y','yCurs_y'};
subj = 2;

clear Hur Hud B F
a = data.rot{subj}.(block);
count = 1;
for j = 1:2
    for k = 1:2
        Hur(:,:,count) = permute(a.Hur(j,k,:,:),[4 3 1 2]);
        Hud(:,:,count) = permute(a.Hud(j,k,:,:),[4 3 1 2]);
        
        B(:,:,count) = permute(a.B(j,k,:,:),[4 3 1 2]);
        F(:,:,count) = permute(a.F(j,k,:,:),[4 3 1 2]);
        count = count + 1;
    end
end

graph_names = {'{xx}','{yx}','{xy}','{yy}'};

figure(1); clf
for j = 1:4
    subplot(2,2,j); hold on
    plot([-1 1],[0 0],'k','HandleVisibility','off')
    plot([0 0],[-1 1],'k','HandleVisibility','off')
%     plot(mean(Hur(:,:,j),1),'-ok','MarkerFaceColor','k')
    for k = 1:Nfreq
        plot(real(Hur(:,k,j)),imag(Hur(:,k,j)),'.','Color',col(k,:),'MarkerSize',20)
    end
    title(['Hur_' graph_names{j}])
    axis([-1.2 1.2 -1.2 1.2])
    axis square
end

figure(2); clf
for j = 1:4
    subplot(2,2,j); hold on
    plot([-1 1],[0 0],'k','HandleVisibility','off')
    plot([0 0],[-1 1],'k','HandleVisibility','off')
%     plot(mean(Hud(:,:,j),1),'-ok','MarkerFaceColor','k')
    for k = 1:Nfreq
        plot(real(Hud(:,k,j)),imag(Hud(:,k,j)),'.','Color',col(k,:),'MarkerSize',20)
    end
    title(['Hud_' graph_names{j}])
    axis([-1.2 1.2 -1.2 1.2])
    axis square
end

figure(3); clf
for j = 1:4
    subplot(2,2,j); hold on
    plot([-1 1],[0 0],'k','HandleVisibility','off')
    plot([0 0],[-1 1],'k','HandleVisibility','off')
%     plot(mean(B(:,:,j),1),'-ok','MarkerFaceColor','k')
    for k = 1:Nfreq
        plot(real(B(:,k,j)),imag(B(:,k,j)),'.','Color',col(k,:),'MarkerSize',20)
    end
    title(['B_' graph_names{j}])
    axis([-1.2 1.2 -1.2 1.2])
    axis square
end

figure(4); clf
for j = 1:4
    subplot(2,2,j); hold on
    plot([-1 1],[0 0],'k','HandleVisibility','off')
    plot([0 0],[-1 1],'k','HandleVisibility','off')
%     plot(mean(F(:,:,j),1),'-ok','MarkerFaceColor','k')
    for k = 1:Nfreq
        plot(real(F(:,k,j)),imag(F(:,k,j)),'.','Color',col(k,:),'MarkerSize',20)
    end
    title(['F_' graph_names{j}])
    axis([-1.2 1.2 -1.2 1.2])
    axis square
end

%%

Nsubj = length(data);
Nfreq = 12;
col = copper;
col = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);
groups = {'rot','mir'};
block = 'B3';

for m = 1:length(groups)
    Nsubj = length(data.(groups{m}));
    
    Hur = NaN(4,Nfreq,4,Nsubj);
    Hud = NaN(4,Nfreq,4,Nsubj);
    B = NaN(4,Nfreq,4,Nsubj);
    F = NaN(4,Nfreq,4,Nsubj);
    
    for i = 1:Nsubj
        a = data.(groups{m}){i}.(block);
        count = 1;
        Ntrials = size(a.F,4);
        for j = 1:2
            for k = 1:2
                Hur(i,:,count,1:Ntrials) = a.Hur(j,k,:,:);
                Hud(i,:,count,1:Ntrials) = a.Hud(j,k,:,:);

                B(i,:,count,1:Ntrials) = a.B(j,k,:,:);
                F(i,:,count,1:Ntrials) = a.F(j,k,:,:);
                count = count + 1;
            end
        end
    end
    
    Hur = mean(Hur,4,'omitnan');
    Hud = mean(Hud,4,'omitnan');
    B = mean(B,4,'omitnan');
    F = mean(F,4,'omitnan');
    
    graph_names = {'{xx}','{yx}','{xy}','{yy}'};
    
    figure((m-1)*4+1); clf
    for j = 1:4
        subplot(2,2,j); hold on
        plot([-1 1],[0 0],'k','HandleVisibility','off')
        plot([0 0],[-1 1],'k','HandleVisibility','off')
        %     plot(mean(Hur(:,:,j),1),'-ok','MarkerFaceColor','k')
        for k = 1:Nfreq
            plot(real(Hur(:,k,j)),imag(Hur(:,k,j)),'.','Color',col(k,:),'MarkerSize',20)
        end
        title(['Hur_' graph_names{j}])
        axis([-1.2 1.2 -1.2 1.2])
        axis square
    end
    
    figure((m-1)*4+2); clf
    for j = 1:4
        subplot(2,2,j); hold on
        plot([-1 1],[0 0],'k','HandleVisibility','off')
        plot([0 0],[-1 1],'k','HandleVisibility','off')
        %     plot(mean(Hud(:,:,j),1),'-ok','MarkerFaceColor','k')
        for k = 1:Nfreq
            plot(real(Hud(:,k,j)),imag(Hud(:,k,j)),'.','Color',col(k,:),'MarkerSize',20)
        end
        title(['Hud_' graph_names{j}])
        axis([-1.2 1.2 -1.2 1.2])
        axis square
    end
    
    figure((m-1)*4+3); clf
    for j = 1:4
        subplot(2,2,j); hold on
        plot([-1 1],[0 0],'k','HandleVisibility','off')
        plot([0 0],[-1 1],'k','HandleVisibility','off')
        %     plot(mean(B(:,:,j),1),'-ok','MarkerFaceColor','k')
        for k = 1:Nfreq
            plot(real(B(:,k,j)),imag(B(:,k,j)),'.','Color',col(k,:),'MarkerSize',20)
        end
        title(['B_' graph_names{j}])
        axis([-1.2 1.2 -1.2 1.2])
        axis square
    end
    
    figure((m-1)*4+4); clf
    for j = 1:4
        subplot(2,2,j); hold on
        plot([-1 1],[0 0],'k','HandleVisibility','off')
        plot([0 0],[-1 1],'k','HandleVisibility','off')
        %     plot(mean(F(:,:,j),1),'-ok','MarkerFaceColor','k')
        for k = 1:Nfreq
            plot(real(F(:,k,j)),imag(F(:,k,j)),'.','Color',col(k,:),'MarkerSize',20)
        end
        title(['F_' graph_names{j}])
        axis([-1.2 1.2 -1.2 1.2])
        axis square
    end
end