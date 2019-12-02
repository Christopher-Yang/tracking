function phasor_1d(data)

subj = 2;

figure(1); clf
subplot(1,2,1); hold on
plot([-2 2],[0 0],'k','LineWidth',1)
plot([0 0],[-2 2],'k','LineWidth',1)
plot(permute(data.org.phasors.cursor.x_x_all.ratio,[2 1]),'.','MarkerSize',20)
axis([-1 1 -1 1])
pbaspect([1 1 1])
title('Original mapping')

subplot(1,2,2); hold on
plot([-2 2],[0 0],'k','LineWidth',1)
plot([0 0],[-2 2],'k','LineWidth',1)
plot(permute(data.new.phasors.cursor.x_x_all.ratio,[2 1]),'.','MarkerSize',20)
axis([-1 1 -1 1])
pbaspect([1 1 1])
title('New mapping')