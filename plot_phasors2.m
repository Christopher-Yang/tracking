subj = 2;

Nfreq = length(data{1}{1}.sineParams.tX_freq{1});
col = copper;
col = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

names1 = {'xTarg_x','yTarg_x','xTarg_y','yTarg_y'};
names2 = {'xCurs_x','yCurs_x','xCurs_y','yCurs_y'};
clear Hur Hud Hur2 Hud2 B B2 F F2
for j = 1:4
    for i = 1:6
        a = data{subj}{i}.phasors;
        Hur(:,i,j) = a.(names1{j}){1}.ratio;
        Hud(:,i,j) = a.(names2{j}){2}.ratio;
        Hur2(:,i,j) = reshape(permute([a.(names1{j}){3}.ratio a.(names1{j}){4}.ratio],[2 1]),[6 1]);
        Hud2(:,i,j) = reshape(permute([a.(names2{j}){4}.ratio a.(names2{j}){3}.ratio],[2 1]),[6 1]);
    end
end

count = 1;
for k = 1:2
    for j = 1:2
        for i = 1:6
            a = data{subj}{i};
            B(:,i,count) = a.B(j,k,:);
            B2(:,i,count) = a.B2(j,k,:);
            F(:,i,count) = a.F(j,k,:);
            F2(:,i,count) = a.F2(j,k,:);
        end
        count = count + 1;
    end
end

clear Hur_trace Hud_trace Hur2_trace Hud2_trace
for i = 1:4
    for j = 1:6
        Hur_trace(j,i) = trace(cov(real(Hur(j,:,i)),imag(Hur(j,:,i))));
        Hud_trace(j,i) = trace(cov(real(Hud(j,:,i)),imag(Hud(j,:,i))));
        Hur2_trace(j,i) = trace(cov(real(Hur2(j,:,i)),imag(Hur2(j,:,i))));
        Hud2_trace(j,i) = trace(cov(real(Hud2(j,:,i)),imag(Hud2(j,:,i))));
    end
end

figure(1); clf
for i = 1:4
    subplot(4,2,(i-1)*2 + 1); hold on
    plot(Hur_trace(:,i),'LineWidth',2)
    plot(Hur2_trace(:,i),'--','LineWidth',2)
    if i == 1
        title('H_{ur}')
    end
    
    subplot(4,2,(i-1)*2 + 2); hold on
    plot(Hud_trace(:,i),'LineWidth',2)
    plot(Hud2_trace(:,i),'--','LineWidth',2)
    if i == 1
        title('H_{ud}')
    end
end

%%

figure(3); clf
for j = 1:4
    subplot(2,2,j); hold on
    plot([-1 1],[0 0],'k','HandleVisibility','off')
    plot([0 0],[-1 1],'k','HandleVisibility','off')
    for i = 1:6
        plot(Hud(i,:,j),'.','MarkerSize',15,'Color',col(i,:))
        plot(nanmean(Hud(i,:,j)),'.','MarkerSize',30,'Color',col(i,:))
    end
    title('Hud')
    axis([-1.3 1.3 -1.3 1.3])
    axis square
end

figure(4); clf
for j = 1:4
    subplot(2,2,j); hold on
    plot([-1 1],[0 0],'k','HandleVisibility','off')
    plot([0 0],[-1 1],'k','HandleVisibility','off')
    for i = 1:6
        plot(Hur(i,:,j),'.','MarkerSize',15,'Color',col(i,:))
        plot(nanmean(Hur(i,:,j)),'.','MarkerSize',30,'Color',col(i,:))
    end
    title('Hur')
    axis([-1.3 1.3 -1.3 1.3])
    axis square
end

figure(5); clf
for j = 1:4
    subplot(2,2,j); hold on
    plot([-1 1],[0 0],'k','HandleVisibility','off')
    plot([0 0],[-1 1],'k','HandleVisibility','off')
    for i = 1:6
        plot(Hud2(i,:,j),'.','MarkerSize',15,'Color',col(i,:))
        plot(nanmean(Hud2(i,:,j)),'.','MarkerSize',30,'Color',col(i,:))
    end
    title('Hud')
    axis([-1.3 1.3 -1.3 1.3])
    axis square
end

figure(6); clf
for j = 1:4
    subplot(2,2,j); hold on
    plot([-1 1],[0 0],'k','HandleVisibility','off')
    plot([0 0],[-1 1],'k','HandleVisibility','off')
    for i = 1:6
        plot(Hur2(i,:,j),'.','MarkerSize',15,'Color',col(i,:))
        plot(nanmean(Hur2(i,:,j)),'.','MarkerSize',30,'Color',col(i,:))
    end
    title('Hur')
    axis([-1.3 1.3 -1.3 1.3])
    axis square
end

%% 

figure(3); clf
for j = 1:4
    subplot(2,2,j); hold on
    plot([-1 1],[0 0],'k','HandleVisibility','off')
    plot([0 0],[-1 1],'k','HandleVisibility','off')
    for i = 1:6
        plot(F(i,:,j),'.','MarkerSize',15,'Color',col(i,:))
        plot(nanmean(F(i,:,j)),'.','MarkerSize',30,'Color',col(i,:))
    end
    title('F')
    axis([-1.3 1.3 -1.3 1.3])
    axis square
end

figure(4); clf
for j = 1:4
    subplot(2,2,j); hold on
    plot([-1 1],[0 0],'k','HandleVisibility','off')
    plot([0 0],[-1 1],'k','HandleVisibility','off')
    for i = 1:6
        plot(B(i,:,j),'.','MarkerSize',15,'Color',col(i,:))
        plot(nanmean(B(i,:,j)),'.','MarkerSize',30,'Color',col(i,:))
    end
    title('B')
    axis([-1.3 1.3 -1.3 1.3])
    axis square
end

figure(5); clf
for j = 1:4
    subplot(2,2,j); hold on
    plot([-1 1],[0 0],'k','HandleVisibility','off')
    plot([0 0],[-1 1],'k','HandleVisibility','off')
    for i = 1:6
        plot(F2(i,:,j),'.','MarkerSize',15,'Color',col(i,:))
        plot(nanmean(F2(i,:,j)),'.','MarkerSize',30,'Color',col(i,:))
    end
    title('F')
    axis([-1.3 1.3 -1.3 1.3])
    axis square
end

figure(6); clf
for j = 1:4
    subplot(2,2,j); hold on
    plot([-1 1],[0 0],'k','HandleVisibility','off')
    plot([0 0],[-1 1],'k','HandleVisibility','off')
    for i = 1:6
        plot(B2(i,:,j),'.','MarkerSize',15,'Color',col(i,:))
        plot(nanmean(B2(i,:,j)),'.','MarkerSize',30,'Color',col(i,:))
    end
    title('B')
    axis([-1.3 1.3 -1.3 1.3])
    axis square
end