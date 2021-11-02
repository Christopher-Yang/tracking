Nsubj = length(data.rot);
Nfreq = 12;
col = copper;
col = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

col1 = [1 0 0];
col2 = [1 1 0];

map = [linspace(col1(1),col2(1),Nfreq)' linspace(col1(2),col2(2),Nfreq)' linspace(col1(3),col2(3),Nfreq)'];

block = 'B3';
names = {'Hur','Hud','B','F'};

clear z
for m = 1:length(names)
    for i = 1:Nsubj
        a = data.rot{i}.(block);
        count = 1;
        for j = 1:2
            for k = 1:2
                z.(names{m})(i,:,count) = mean(a.(names{m})(j,k,:,:),4);
                count = count + 1;
            end
        end
    end
end

graph_names = {'{xx}','{yx}','{xy}','{yy}'};

for i = 1:length(names)
    figure(i); clf
    for j = 1:4
        subplot(2,2,j); hold on
        plot([-1 1],[0 0],'k','HandleVisibility','off')
        plot([0 0],[-1 1],'k','HandleVisibility','off')
        for k = 1:Nfreq
            plot(real(z.(names{i})(:,k,j)),imag(z.(names{i})(:,k,j)),'.','Color',map(k,:),'MarkerSize',10)
        end
        plot(mean(z.(names{i})(:,:,j),1),'k','LineWidth',1)
        for k = 1:Nfreq
            plot(mean(z.(names{i})(:,k,j),1),'-ko','LineWidth',1,'MarkerFaceColor',map(k,:))
        end
        title([names{i} '_' graph_names{j}])
        axis([-1.2 1.2 -1.2 1.2])
        axis square
    end
end

%%
figure(5); clf
for j = 1:4
    subplot(2,2,j); hold on
    plot(abs(z.F(:,:,j))','Color',[0.3 0.3 0.3])
    plot(abs(mean(z.F(:,:,j),1)),'k','LineWidth',2)
end
