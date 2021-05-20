k = NaN(7,10);

x = 3; % axis
y = 5; % block
q = 2; % group

for p = 1:10
    k(:,p) = mean(squeeze(fracOpt(x,y,:,:,p,q)),2);
end

mean(k,2)

%%
figure(4); clf; hold on
histogram(fracOpt([2 3],5,:,:,:,:),0:0.1:1,'Normalization','probability')
axis([0 1 0 1])
ylabel('Probability')
xlabel('Proportion of gain retained')

sum(sum(sum(sum(sum(sum(fracOpt([2 3],5,:,:,:,:)>0.5))))))/numel(fracOpt([2 3],5,:,:,:,:));

%%
x = 3; % axis
y = 5; % block
m = 1; % frequency
p = 1; % subject
q = 2; % group

phasor1 = reshape(thetaOpt([1 2],y,1,:,p,q),[12 1]);

figure(7); clf; hold on
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
% plot

%%
for k = 1:2
    if k == 1
        template = data.(groups{q}){p}.baseline.Rhand.phasors;
    else
        template = data.(groups{q}){p}.pert4.Rhand.phasors;
    end
    % set variables for  plotting
    Nfreq = 7;
    
    % set color map
    col1 = [1 0.9 0.3];
    col2 = [1 0 0];
    colors = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2)...
        ,Nfreq)', linspace(col1(3),col2(3),Nfreq)'];
    
    % draw x- and y-axes and label plots
    figure(4+k); clf; hold on
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
        plot(template.x_x.ratio(i,:),'.','Color',colors(i,:),'MarkerSize',5)
        error_ellipse(cov(real(template.x_x.ratio(i,:)),imag(...
            template.x_x.ratio(i,:))),[real(mean(template.x_x.ratio(i,:))); ...
            imag(mean(template.x_x.ratio(i,:)))],'conf',0.95,'color'...
            ,colors(i,:));
        if(i>1)
            plot(mean(template.x_x.ratio(i-1:i,:),2),'.-','Color',...
                colors(i-1,:),'MarkerSize',10,'linewidth',2)
        end
        
        subplot(2,2,2)
        plot(template.y_x.ratio(i,:),'.','Color',colors(i,:),'MarkerSize',5)
        error_ellipse(cov(real(template.y_x.ratio(i,:)),imag(...
            template.y_x.ratio(i,:))),[real(mean(template.y_x.ratio(i,:))); ...
            imag(mean(template.y_x.ratio(i,:)))],'conf',0.95,'color'...
            ,colors(i,:));
        if(i>1)
            plot(mean(template.y_x.ratio(i-1:i,:),2),'.-','Color'...
                ,colors(i-1,:),'MarkerSize',10,'linewidth',2)
        end
        
        subplot(2,2,3)
        plot(template.x_y.ratio(i,:),'.','Color',colors(i,:),'MarkerSize',5)
        error_ellipse(cov(real(template.x_y.ratio(i,:)),imag(...
            template.x_y.ratio(i,:))),[real(mean(template.x_y.ratio(i,:))); ...
            imag(mean(template.x_y.ratio(i,:)))],'conf',0.95,'color'...
            ,colors(i,:));
        if(i>1)
            plot(mean(template.x_y.ratio(i-1:i,:),2),'.-','Color'...
                ,colors(i-1,:),'MarkerSize',10,'linewidth',2)
        end
        
        subplot(2,2,4)
        plot(template.y_y.ratio(i,:),'.','Color',colors(i,:),'MarkerSize',5)
        error_ellipse(cov(real(template.y_y.ratio(i,:)),imag(...
            template.y_y.ratio(i,:))),[real(mean(template.y_y.ratio(i,:))); ...
            imag(mean(template.y_y.ratio(i,:)))],'conf',0.95,'color'...
            ,colors(i,:));
        if(i>1)
            plot(mean(template.y_y.ratio(i-1:i,:),2),'.-','Color'...
                ,colors(i-1,:),'MarkerSize',10,'linewidth',2)
        end
    end
end