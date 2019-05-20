function graph_phasor(data, subj)

Nfreq = 7;
% col1 = [1 0.83 0.33];
% col2 = [1 0 0];
% colors = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2),Nfreq)', linspace(col1(3),col2(3),Nfreq)'];

rat = data.L.(subj).B1.x_x_all.ratio;
rat(:,:,2) = data.R.(subj).B1.x_x_all.ratio;

a = mean(rat,2);
mu(:,1,:) = real(a);
mu(:,2,:) = imag(a);

figure
subplot(1,2,1); hold on
plot([-1.5 1.5],[0 0],'k')
plot([0 0],[-1.5 1.5],'k')
for i = 1:Nfreq
    plot(rat(i,:,1),'.')
    error_ellipse(cov(real(rat(i,:,1)),imag(rat(i,:,1))),mu(i,:,1),'conf',0.95,'style','k');
end
pbaspect([1 1 1])
title('Left Hand')

subplot(1,2,2); hold on
plot([-1.5 1.5],[0 0],'k')
plot([0 0],[-1.5 1.5],'k')
for i = 1:Nfreq
    plot(rat(i,:,2),'.')
    error_ellipse(cov(real(rat(i,:,1)),imag(rat(i,:,1))),mu(i,:,2),'conf',0.95,'style','k');
end
pbaspect([1 1 1])
title('Right Hand')

end