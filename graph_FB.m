subj = 1;
block = 1;
% phasor = reshape(permute(data{subj}{block}.B,[2 1 3]),[4 6]);
Nfreq = length(data{1}{1}.sineParams.tX_freq{1});
col = copper;
col = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

figure(1); clf
phasor = reshape(permute(data{subj}{block}.B,[2 1 3]),[4 6]);
phasor2 = reshape(permute(data{subj}{block}.B2,[2 1 3]),[4 6]);
for i = 1:4
    subplot(2,2,i); hold on
    plot([-1 1],[0 0],'k','HandleVisibility','off')
    plot([0 0],[-1 1],'k','HandleVisibility','off')
    for j = 1:6
        plot(phasor(i,j),'.','Color',col(j,:),'MarkerSize',20)
        plot(phasor2(i,j),'.','Color',col(j,:),'MarkerSize',20)
    end
    plot(phasor(i,:),'k')
    plot(phasor2(i,:),'--k')
    axis equal
end

%%
Nfreq = length(data{1}{1}.sineParams.tX_freq{1});
col = copper;
col = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

figure(2); clf
for j = [1 2 6]
    gain = abs(reshape(permute(data{subj}{j}.B2,[2 1 3]),[4 6]));
    if j == 1
        gblocks = 1:4;
    else
        gblocks = [2 1 4 3];
    end
    for i = 1:4
        subplot(2,2,i); hold on
        plot(gain(gblocks(i),:),'LineWidth',2)
        axis([1 Nfreq -0.1 2])
    end
end
legend({'Baseline','Early','Late'},'Location','northwest')

%%
figure(3); clf
for j = [2 6]
    phase = angle(reshape(permute(data{subj}{j}.B,[2 1 3]),[4 6]));
    phase = unwrap(phase,[],2)*180/pi;
    for i = 1:4
        subplot(2,2,i); hold on
        plot(phase(i,:),'Color',col(j,:))
%         axis([1 Nfreq -0.1 1.5])
    end
end