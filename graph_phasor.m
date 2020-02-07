function graph_phasor(data,f)
% plots the phasors for a given subject and a given block, specified by the
% input, data

if(nargin<2)
    f = 1; % default plot as figure 1
end

% set variables for  plotting
Nfreq = 7;

% set color map
col1 = [1 0.9 0.3];
col2 = [1 0 0];
colors = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2)...
    ,Nfreq)', linspace(col1(3),col2(3),Nfreq)'];

% draw x- and y-axes and label plots
figure(f); clf; hold on
subplot(2,2,1); hold on
plot([-1.5 1.5],[0 0],'k')
plot([0 0],[-1.5 1.5],'k')
axis square
axis([-1.25 1.25 -1.25 1.25])
title('X_{target} -> X_{hand}')
xticks(-2:1:2)
yticks(-2:1:2)

subplot(2,2,2); hold on
plot([-1.5 1.5],[0 0],'k')
plot([0 0],[-1.5 1.5],'k')
axis square
axis([-1.25 1.25 -1.25 1.25])
title('Y_{target} -> X_{hand}')
xticks(-2:1:2)
yticks(-2:1:2)

subplot(2,2,3); hold on
plot([-1.5 1.5],[0 0],'k')
plot([0 0],[-1.5 1.5],'k')
axis square
axis([-1.25 1.25 -1.25 1.25])
title('X_{target} -> Y_{hand}')
xticks(-2:1:2)
yticks(-2:1:2)

subplot(2,2,4); hold on
plot([-1.5 1.5],[0 0],'k')
plot([0 0],[-1.5 1.5],'k')
axis square
axis([-1.25 1.25 -1.25 1.25])
title('Y_{target} -> Y_{hand}')
xticks(-2:1:2)
yticks(-2:1:2)

% plot phasors
for i = 1:Nfreq
    subplot(2,2,1)
    plot(data.x_x_all.ratio(i,:),'.','Color',colors(i,:),'MarkerSize',5)
    error_ellipse(cov(real(data.x_x_all.ratio(i,:)),imag(...
        data.x_x_all.ratio(i,:))),[real(mean(data.x_x_all.ratio(i,:))); ...
        imag(mean(data.x_x_all.ratio(i,:)))],'conf',0.95,'color'...
        ,colors(i,:));
    if(i>1)
        plot(mean(data.x_x_all.ratio(i-1:i,:),2),'.-','Color',...
            colors(i-1,:),'MarkerSize',10,'linewidth',2)
    end
    
    subplot(2,2,2)
    plot(data.y_x_all.ratio(i,:),'.','Color',colors(i,:),'MarkerSize',5)
    error_ellipse(cov(real(data.y_x_all.ratio(i,:)),imag(...
        data.y_x_all.ratio(i,:))),[real(mean(data.y_x_all.ratio(i,:))); ...
        imag(mean(data.y_x_all.ratio(i,:)))],'conf',0.95,'color'...
        ,colors(i,:));
    if(i>1)
        plot(mean(data.y_x_all.ratio(i-1:i,:),2),'.-','Color'...
            ,colors(i-1,:),'MarkerSize',10,'linewidth',2)
    end
    
    subplot(2,2,3)
    plot(data.x_y_all.ratio(i,:),'.','Color',colors(i,:),'MarkerSize',5)
    error_ellipse(cov(real(data.x_y_all.ratio(i,:)),imag(...
        data.x_y_all.ratio(i,:))),[real(mean(data.x_y_all.ratio(i,:))); ...
        imag(mean(data.x_y_all.ratio(i,:)))],'conf',0.95,'color'...
        ,colors(i,:));
    if(i>1)
        plot(mean(data.x_y_all.ratio(i-1:i,:),2),'.-','Color'...
            ,colors(i-1,:),'MarkerSize',10,'linewidth',2)
    end
    
    subplot(2,2,4)
    plot(data.y_y_all.ratio(i,:),'.','Color',colors(i,:),'MarkerSize',5)
    error_ellipse(cov(real(data.y_y_all.ratio(i,:)),imag(...
        data.y_y_all.ratio(i,:))),[real(mean(data.y_y_all.ratio(i,:))); ...
        imag(mean(data.y_y_all.ratio(i,:)))],'conf',0.95,'color'...
        ,colors(i,:));
    if(i>1)
        plot(mean(data.y_y_all.ratio(i-1:i,:),2),'.-','Color'...
            ,colors(i-1,:),'MarkerSize',10,'linewidth',2)
    end
end

set(gcf,'Renderer','painters')