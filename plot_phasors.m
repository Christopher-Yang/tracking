Nsubj = length(data);
Nfreq = 12;
col = copper;
col = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

block = 5;
names1 = {'xTarg_x','yTarg_x','xTarg_y','yTarg_y'};
names2 = {'xCurs_x','yCurs_x','xCurs_y','yCurs_y'};
subj = 1;

clear Hur Hur2 Hud Hud2 B B2 F F2
a = data{subj}{block};
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

graph_names = {'X-->X','Y-->X','X-->Y','Y-->Y'};

figure(1); clf
for j = 1:4
    subplot(2,2,j); hold on
    plot([-1 1],[0 0],'k','HandleVisibility','off')
    plot([0 0],[-1 1],'k','HandleVisibility','off')
%     plot(mean(Hur(:,:,j),1),'-ok','MarkerFaceColor','k')
    for k = 1:Nfreq
        plot(real(Hur(:,k,j)),imag(Hur(:,k,j)),'.','Color',col(k,:),'MarkerSize',20)
    end
    title(graph_names{j})
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
    title(graph_names{j})
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
    title(graph_names{j})
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
    title(graph_names{j})
    axis([-1.2 1.2 -1.2 1.2])
    axis square
end
