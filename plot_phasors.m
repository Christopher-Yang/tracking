Nsubj = length(data);
Nfreq = length(data{1}{1}.sineParams.tX_freq{1});
col = copper;
col = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

block = 1;
names1 = {'xTarg_x','yTarg_x','xTarg_y','yTarg_y'};
names2 = {'xCurs_x','yCurs_x','xCurs_y','yCurs_y'};

for i = 1:Nsubj
    a = data{i}{block}.phasors;
    for j = 1:4
        phasors1{j}(i,:) = a.(names1{j}){1}.ratio;
        phasors2{j}(i,:) = a.(names2{j}){2}.ratio;
        phasors3{j}(i,:) = reshape(permute([a.(names1{j}){3}.ratio a.(names1{j}){4}.ratio],[2 1]),[1 6]);
        phasors4{j}(i,:) = reshape(permute([a.(names2{j}){4}.ratio a.(names2{j}){3}.ratio],[2 1]),[1 6]);
    end
end

graph_names = {'X-->X','Y-->X','X-->Y','Y-->Y'};

figure(1); clf
for j = 1:4
    subplot(2,2,j); hold on
    plot([-1 1],[0 0],'k','HandleVisibility','off')
    plot([0 0],[-1 1],'k','HandleVisibility','off')
    plot(phasors1{j},'k')
    plot(phasors3{j},'--k')
    for k = 1:Nfreq
        plot(phasors1{j}(:,k),'.','Color',col(k,:),'MarkerSize',20)
        plot(phasors3{j}(:,k),'.','Color',col(k,:),'MarkerSize',20)
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
    plot(phasors2{j},'k')
    plot(phasors4{j},'--k')
    for k = 1:Nfreq
        plot(phasors2{j}(:,k),'.','Color',col(k,:),'MarkerSize',20)
        plot(phasors4{j}(:,k),'.','Color',col(k,:),'MarkerSize',20)
    end
    title(graph_names{j})
    if j == 1
        legend({'Individual','Dual'})
    end
    axis([-1.2 1.2 -1.2 1.2])
    axis square
end
