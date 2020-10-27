Nsubj = length(data);
Nfreq = length(data{1}{1}.sineParams.tX_freq{1});
col = copper;
col = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

block = 1;
names1 = {'xTarg_x','yTarg_x','xTarg_y','yTarg_y'};
names2 = {'xCurs_x','yCurs_x','xCurs_y','yCurs_y'};
subjs = 1;

clear Hur Hur2 Hud Hud2 B B2 F F2
for i = subjs
    a = data{i}{block};
    for j = 1:4
        Hur(i,:,j) = a.phasors.(names1{j}){1}.ratio;
        Hud(i,:,j) = a.phasors.(names2{j}){2}.ratio;
%         Hur2(i,:,j) = reshape(permute([a.phasors.(names1{j}){3}.ratio a.phasors.(names1{j}){4}.ratio],[2 1]),[1 6]);
%         Hud2(i,:,j) = reshape(permute([a.phasors.(names2{j}){4}.ratio a.phasors.(names2{j}){3}.ratio],[2 1]),[1 6]);
    end
    
    count = 1;
    for j = 1:2
        for k = 1:2
            B(i,:,count) = a.B(j,k,:);
            B2(i,:,count) = a.B2(j,k,:);
            F(i,:,count) = a.F(j,k,:);
            F2(i,:,count) = a.F2(j,k,:);
            count = count + 1;
        end
    end
end

graph_names = {'X-->X','Y-->X','X-->Y','Y-->Y'};

figure(1); clf
for j = 1:4
    subplot(2,2,j); hold on
    plot([-1 1],[0 0],'k','HandleVisibility','off')
    plot([0 0],[-1 1],'k','HandleVisibility','off')
    plot(mean(Hur(:,:,j),1),'-ok','MarkerFaceColor','k')
    plot(mean(Hur2(:,:,j),1),'--ok','MarkerFaceColor','k')
    for k = 1:Nfreq
        plot(Hur(:,k,j),'.','Color',col(k,:),'MarkerSize',20)
        plot(Hur2(:,k,j),'.','Color',col(k,:),'MarkerSize',20)
    end
    title(graph_names{j})
    if j == 1
        legend({'Individual','Dual'})
    end
    axis([-1.2 1.2 -1.2 1.2])
    axis square
end

figure(2); clf
for j = 1:4
    subplot(2,2,j); hold on
    plot([-1 1],[0 0],'k','HandleVisibility','off')
    plot([0 0],[-1 1],'k','HandleVisibility','off')
    plot(mean(Hud(:,:,j),1),'-ok','MarkerFaceColor','k')
    plot(mean(Hud2(:,:,j),1),'--ok','MarkerFaceColor','k')
    for k = 1:Nfreq
        plot(Hud(:,k,j),'.','Color',col(k,:),'MarkerSize',20)
        plot(Hud2(:,k,j),'.','Color',col(k,:),'MarkerSize',20)
    end
    title(graph_names{j})
    if j == 1
        legend({'Individual','Dual'})
    end
    axis([-1.2 1.2 -1.2 1.2])
    axis square
end

% figure(3); clf
% for j = 1:4
%     subplot(2,2,j); hold on
%     plot([-1 1],[0 0],'k','HandleVisibility','off')
%     plot([0 0],[-1 1],'k','HandleVisibility','off')
%     plot(mean(B(:,:,j),1),'-ok','MarkerFaceColor','k')
%     plot(mean(B2(:,:,j),1),'--ok','MarkerFaceColor','k')
%     for k = 1:Nfreq
%         plot(B(:,k,j),'.','Color',col(k,:),'MarkerSize',20)
%         plot(B2(:,k,j),'.','Color',col(k,:),'MarkerSize',20)
%     end
%     title(graph_names{j})
%     if j == 1
%         legend({'Individual','Dual'})
%     end
%     axis([-1.2 1.2 -1.2 1.2])
%     axis square
% end
% 
% figure(4); clf
% for j = 1:4
%     subplot(2,2,j); hold on
%     plot([-1 1],[0 0],'k','HandleVisibility','off')
%     plot([0 0],[-1 1],'k','HandleVisibility','off')
%     plot(mean(F(:,:,j),1),'-ok','MarkerFaceColor','k')
%     plot(mean(F2(:,:,j),1),'--ok','MarkerFaceColor','k')
%     for k = 1:Nfreq
%         plot(F(:,k,j),'.','Color',col(k,:),'MarkerSize',20)
%         plot(F2(:,k,j),'.','Color',col(k,:),'MarkerSize',20)
%     end
%     title(graph_names{j})
%     if j == 1
%         legend({'Individual','Dual'})
%     end
%     axis([-1.2 1.2 -1.2 1.2])
%     axis square
% end
