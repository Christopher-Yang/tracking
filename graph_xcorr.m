function graph_xcorr(data,block_name,graph_name)

col = lines;
col = col(1:7,:);
idx = find(contains(graph_name,'('));

Nsubj = length(data)-1;
Nblock = length(fieldnames(data{1}));

for i = 1:Nsubj
    for j = 1:Nblock
        xCor_all(:,j,i) = data{i}.(block_name{j}).xCor;
        yCor_all(:,j,i) = data{i}.(block_name{j}).yCor;
    end
end

xCor = mean(xCor_all,3);
yCor = mean(yCor_all,3);

n = 1:numel(xCor);
n = reshape(n,[size(xCor,1) Nblock]);
n2 = setdiff(1:Nblock,idx);

figure(1); clf; hold on
plot([-10 -9],[0 0],'LineWidth',1.5)
plot([-10 -9],[0 0],'LineWidth',1.5)
plot([0 n(end)],[0 0],'--k','LineWidth',1)
plot(n(:,n2),xCor(:,n2),'Color',col(1,:),'LineWidth',1.5)
plot(n(:,idx),xCor(:,idx),'-o','Color',col(1,:),'LineWidth',2)
plot(n(:,n2),yCor(:,n2),'Color',col(2,:),'LineWidth',1.5)
plot(n(:,idx),yCor(:,idx),'-o','Color',col(2,:),'LineWidth',2)
for i = 1:Nsubj
    plot(n,xCor_all(:,:,i),'Color',[col(1,:) 0.4])
    plot(n,yCor_all(:,:,i),'Color',[col(2,:) 0.4])
end
xlim([1 n(end)])
legend({'X','Y'});
xticks(idx*size(xCor,1)-2)
xticklabels(graph_name(idx))

end