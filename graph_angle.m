function graph_angle(data, groups, gblocks, theta)

    col = lines;
    col = col(1:7,:);
    group_names = {'Rotation','Mirror Reversal'};
    freqs_x = {'0.10 Hz','0.25 Hz','0.55 Hz','0.85 Hz','1.15 Hz','1.55 Hz','2.05 Hz'};
    freqs_y = {'0.15 Hz','0.35 Hz','0.65 Hz','0.95 Hz','1.45 Hz','1.85 Hz','2.15 Hz'};
    f_x = data.(groups{1}).avg.x_x.freqs;
    f_y = data.(groups{1}).avg.y_y.freqs;
    leg = {'Low frequencies','High frequencies'};
    line_width = 1.25;
    font_size = 12;
    Nblock = size(data.(groups{1}).avg.x_x.fft,1);
    theta_rad = theta*pi/180;
    
%     a = figure('DefaultTextFontName','Arial','DefaultAxesFontName','Arial');
    b = figure('Name','Angle of Movement Error','NumberTitle','off','DefaultTextFontName','Arial','DefaultAxesFontName','Arial');
%     c = figure('DefaultTextFontName','Arial','DefaultAxesFontName','Arial');

    for i = 1:length(groups)
        dat = data.(groups{i}).avg;
        boot = NaN(length(f_x),1000);
        bootX = NaN(Nblock,length(f_x),1000);
        bootY = NaN(Nblock,length(f_y),1000);
        err = NaN(length(f_x),2);
        errX = NaN(Nblock,length(f_x),2);
        errY = NaN(Nblock,length(f_y),2);
        if strcmp(groups{i},'rot')
            angX = atan2(mean(dat.x_y.all_amp,3),mean(dat.x_x.all_amp,3));
            angY = atan2(mean(dat.y_x.all_amp,3),mean(dat.y_y.all_amp,3));
            angX_all = atan2(dat.x_y.all_amp, dat.x_x.all_amp);
            angY_all = atan2(dat.y_x.all_amp, dat.y_y.all_amp);
        elseif strcmp(groups{i},'rot_i')
            angX = atan2(mean(dat.x_y.all_amp,3),mean(dat.x_x.all_amp,3));
            angY = atan2(mean(dat.y_x.all_amp,3),mean(dat.y_y.all_amp,3));
            angX_all = atan2(dat.x_y.all_amp, dat.x_x.all_amp);
            angY_all = atan2(dat.y_x.all_amp, dat.y_y.all_amp);
        end
        
        angX(2:5,:) = theta_rad - angX(2:5,:);
        angY(2:5,:) = theta_rad - angY(2:5,:);
        angX_all(2:5,:,:) = theta_rad - angX_all(2:5,:,:);
        angY_all(2:5,:,:) = theta_rad - angY_all(2:5,:,:);
        d(:,:,:,1,i) = angX_all([1 6],:,:);
        d(:,:,:,2,i) = angY_all([1 6],:,:);
        
        totalX_all = squeeze(angX_all(5,:,:) - angX_all(1,:,:));
        totalY_all = squeeze(angY_all(5,:,:) - angY_all(1,:,:));
        afterX_all = squeeze(angX_all(6,:,:) - angX_all(1,:,:));
        afterY_all = squeeze(angY_all(6,:,:) - angY_all(1,:,:));
        aftereffX = filloutliers(afterX_all./totalX_all,NaN,2);
        aftereffY = filloutliers(afterY_all./totalY_all,NaN,2);
        aftereff = nanmean(cat(3,aftereffX,aftereffY),3);
        avg = nanmean(aftereff,2);
        
        for j = 1:1000
            p = datasample(angX_all,10,3);
            q = datasample(angY_all,10,3);
            bootX(:,:,j) = mean(p,3);
            bootY(:,:,j) = mean(q,3);
            
            p = datasample(aftereff,10,2);
            boot(:,j) = nanmean(p,2);
        end
        
        boot = sort(boot,2);
        err_p = cat(2,boot(:,26),boot(:,end-25));
        err(:,1) = avg - err_p(:,1);
        err(:,2) = err_p(:,2) - avg;
        err = err';
        err = [err(2,:); err(1,:)];
        
        bootX = sort(bootX,3);
        bootY = sort(bootY,3);
        errX_p = cat(3,bootX(:,:,26),bootX(:,:,end-25));
        errY_p = cat(3,bootY(:,:,26),bootY(:,:,end-25));
        errX(:,:,1) = angX - errX_p(:,:,1);
        errX(:,:,2) = errX_p(:,:,2) - angX;
        errY(:,:,1) = angY - errY_p(:,:,1);
        errY(:,:,2) = errY_p(:,:,2) - angY;
        errX = permute(errX*2/pi,[3 2 1]);
        errY = permute(errY*2/pi,[3 2 1]);
        errX = [errX(2,:,:); errX(1,:,:)];
        errY = [errY(2,:,:); errY(1,:,:)];
        
%% Aftereffects summary
%         figure(a);
%         col2 = [0 0.5 0; 0 0.75 0.75];
%         s = shadedErrorBar(1:7,100*avg,100*err,'lineProps','-o'); hold on;
%         editErrorBar(s,col2(i,:),3);
%         plot(-5:10,zeros(16,1),'-k','HandleVisibility','off');
%         set(gca,'FontSize',16,'TickDir','out','box','off','LineWidth',1,'ytick',-100:25:100,'Xtick',[1 7],'xticklabel',{'Lowest frequency','Highest frequency'})
%         title('Aftereffects','FontSize',font_size);
%         ylabel('Aftereffects (%)','FontSize',font_size);
%         axis([1 7 -30 75]);
%         if i == 2
%             legend(group_names,'Location','northwest','FontSize',18);
%         end
%         legend boxoff;
%% Compensation angle not flipped
        figure(b)
        subplot(2,2,i)
        for j = 1:length(gblocks)
            s = shadedErrorBar(f_x,angX(gblocks(j),:)*2/pi,errX(:,:,gblocks(j)),'lineProps','-o'); hold on; % using 2/pi because 2/pi = 180/pi * 1/90
            editErrorBar(s,col(j,:),line_width);
        end
        plot([0.001 10],[1 1],'--k','LineWidth',2)
        set(gca,'Xscale','log','Ytick',0:1,'FontSize',14,'Xgrid','on','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','TickDir','out','LineWidth',1)
        title(group_names{i},'FontSize',font_size)
        if i == 1
            ylabel('X Compensation Index','FontSize',font_size)
        else
            legend({'Baseline', 'Naive', 'Max training', 'Aftereffects','Ideal'},'Position',[0.80 0.76 0.1 0.1],'FontSize',font_size)
        end
        axis([0.09 2.3 0 1.05])
        pbaspect([1 1 1])
        
        subplot(2,2,i+2)
        for j = 1:length(gblocks)
            s = shadedErrorBar(f_y,angY(gblocks(j),:)*2/pi,errY(:,:,gblocks(j)),'lineProps','-o'); hold on;
            editErrorBar(s,col(j,:),line_width);
        end
        plot([0.001 10],[1 1],'--k','LineWidth',2)
        set(gca,'Xscale','log','Ytick',0:1,'FontSize',14,'Xgrid','on','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','TickDir','out','LineWidth',1)
        xlabel('Frequency (Hz)','FontSize',font_size)
        if i == 1
            ylabel('Y Compensation Index','FontSize',font_size)
        end
        axis([0.09 2.3 0 1.05])
        pbaspect([1 1 1])
%% Compensation angle flipped
%         figure(c)
%         subplot(2,2,i);
%         s = shadedErrorBar(f_x,angX(gblocks(1),:)*180/pi,errX(:,:,gblocks(1)),'lineProps','-o'); hold on;
%         editErrorBar(s,col(1,:),line_width);
%         s = shadedErrorBar(f_y,angY(gblocks(2),:)*180/pi,errY(:,:,gblocks(2)),'lineProps','-o');
%         editErrorBar(s,col(2,:),line_width);
%         s = shadedErrorBar(f_y,angY(gblocks(3),:)*180/pi,errY(:,:,gblocks(3)),'lineProps','-o');
%         editErrorBar(s,col(3,:),line_width);
%         s = shadedErrorBar(f_x,angX(gblocks(4),:)*180/pi,errX(:,:,gblocks(4)),'lineProps','-o');
%         editErrorBar(s,col(4,:),line_width);
%         plot([0.001 10],[90 90],'--k','LineWidth',line_width,'HandleVisibility','off');
%         set(gca,'Xscale','log','Ytick',0:30:90,'FontSize',12,'Xgrid','on','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','TickDir','out','LineWidth',1);
%         title(group_names{i},'FontSize',18);
%         if i == 1
%             ylabel(['X compensation angle (',char(176),')'],'FontSize',12);
%             legend({'Baseline', 'Naive', 'Max training', 'Aftereffects'},'Position',[0.14 0.78 0.1 0.1],'FontSize',12); 
%         end
%         axis([0.09 2.3 0 95]);
%         
%         subplot(2,2,i+2);
%         s = shadedErrorBar(f_y,angY(gblocks(1),:)*180/pi,errY(:,:,gblocks(1)),'lineProps','-o'); hold on;
%         editErrorBar(s,col(1,:),line_width);
%         s = shadedErrorBar(f_x,angX(gblocks(2),:)*180/pi,errX(:,:,gblocks(2)),'lineProps','-o');
%         editErrorBar(s,col(2,:),line_width);
%         s = shadedErrorBar(f_x,angX(gblocks(3),:)*180/pi,errX(:,:,gblocks(3)),'lineProps','-o');
%         editErrorBar(s,col(3,:),line_width);
%         s = shadedErrorBar(f_y,angY(gblocks(4),:)*180/pi,errY(:,:,gblocks(4)),'lineProps','-o');
%         editErrorBar(s,col(4,:),line_width);
%         plot([0.001 10],[90 90],'--k','LineWidth',line_width);
%         set(gca,'Xscale','log','Ytick',0:30:90,'FontSize',12,'Xgrid','on','GridAlpha',0.4,'MinorGridAlpha',0.5,'box','off','TickDir','out','LineWidth',1);
%         xlabel('Frequency (Hz)','FontSize',12);
%         if i == 1
%             ylabel(['Y compensation angle (',char(176),')'],'FontSize',12);
%         end
%         axis([0.09 2.3 0 95]);
    end
end