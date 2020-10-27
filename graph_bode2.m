% plot Bode plots for learning on myself

freqX = data{1}{1}.sineParams.tX_freq{1};
freqY = data{1}{1}.sineParams.tY_freq{1};
Nblock = 6;

Nfreq = length(data{1}{1}.sineParams.tX_freq{1});
col = copper;
col = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

clear B B2 F F2
for i = 1:Nblock
    a = data{1}{i};
    count = 1;
    for j = 1:2
        for k = 1:2
            B(:,i,count) = a.B(j,k,:);
            B2(:,i,count) = a.B2(j,k,:);
            F(:,i,count) = a.F(j,k,:);
            F2(:,i,count) = a.F2(j,k,:);
            count = count + 1;
        end
    end
end

freq = 4;

figure(1); clf
for i = 1:4
    subplot(2,2,i); hold on
    plot([-1 1],[0 0],'k')
    plot([0 0],[-1 1],'k')
    for j = 1:6
        plot(permute(B(freq,j,i),[2 1]),'.','MarkerSize',20,'Color',col(j,:))
    end
end

figure(2); clf
for i = 1:4
    subplot(2,2,i); hold on
    plot([-1 1],[0 0],'k')
    plot([0 0],[-1 1],'k')
    for j = 1:6
        plot(permute(F(freq,j,i),[2 1]),'.','MarkerSize',20,'Color',col(j,:))
    end
end

%%
B_gain = abs(B);
F_gain = abs(F);
B_phase = unwrap(angle(B));
F_phase = unwrap(angle(F));
idx = 1:6;

B_phase(:,5,2) = B_phase(:,5,2) + 2*pi;
B_phase(5:6,6,2) = B_phase(5:6,6,2) - 2*pi;

figure(1); clf
subplot(2,2,1); hold on
for i = 1:Nblock
    plot(freqX(idx),B_gain(idx,i,2),'Color',col(i,:),'LineWidth',2)
end
title('B')
ylabel('Gain')

subplot(2,2,2); hold on
for i = 1:Nblock
    plot(freqX(idx),F_gain(idx,i,2),'Color',col(i,:),'LineWidth',2)
end
title('F')

subplot(2,2,3)
for i = 1:Nblock; hold on
    plot(freqX,B_phase(:,i,2)*180/pi,'Color',col(i,:),'LineWidth',2)
end
ylabel('Phase')

subplot(2,2,4); hold on
for i = 1:Nblock
    plot(freqX,F_phase(:,i,2)*180/pi,'Color',col(i,:),'LineWidth',2)
end

% figure(1); clf; hold on
% plot([-1 1],[0 0],'k','HandleVisibility','off')
% plot([0 0],[-1 1],'k','HandleVisibility','off')
% plot(B(:,1,1))
