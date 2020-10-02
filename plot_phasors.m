Nsubj = length(data);
Nfreq = length(data{1}.sineParams.tX_freq{1});
col = copper;
col = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

for i = 1:Nsubj
    a = data{i}.phasors;
    phasors1(i,:) = a.xTarg_x{1}.ratio;
    phasors2(i,:) = a.xCurs_x{3}.ratio;
    phasors3(i,:) = reshape(permute([a.xTarg_x{5}.ratio a.xTarg_x{7}.ratio],[2 1]),[1 6]);
    phasors4(i,:) = reshape(permute([a.xCurs_x{7}.ratio a.xCurs_x{5}.ratio],[2 1]),[1 6]);
end

subj = 1;
a = data{subj}.phasors;

figure(1); clf
subplot(2,2,1); hold on
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
for k = 1:Nfreq
    plot(phasors1(:,k),'.','Color',col(k,:),'MarkerSize',20)
end
title('Target response')
ylabel('Single perturbation')
axis([-1.2 1.2 -1.2 1.2])
axis square

subplot(2,2,2); hold on
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
for k = 1:Nfreq
    plot(phasors2(:,k),'.','Color',col(k,:),'MarkerSize',20)
end
title('Cursor response')
axis([-1.2 1.2 -1.2 1.2])
axis square

subplot(2,2,3); hold on
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
for k = 1:Nfreq
    plot(phasors3(:,k),'.','Color',col(k,:),'MarkerSize',20)
end
ylabel('Dual perturbation')
axis([-1.2 1.2 -1.2 1.2])
axis square

subplot(2,2,4); hold on
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
for k = 1:Nfreq
    plot(phasors4(:,k),'.','Color',col(k,:),'MarkerSize',20)
end
axis([-1.2 1.2 -1.2 1.2])
axis square