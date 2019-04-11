function graph_bode(data, graph_name, groups)
    
    rng(1);
    titles = {'rotation', 'mirror reversal'};
    Nblock = size(data.rot.avg.x_x.fft,1);
    Nfreq = size(data.rot.avg.x_x.freqs,1);
    
    a = reshape(datasample([1 -1],Nblock*Nfreq),[Nfreq Nblock]);
    x = rand(Nfreq,Nblock).*a*0.01;
    y = rand(Nfreq,Nblock).*a*0.01;
    scale = repmat(1:Nfreq,[Nblock,1])';
    x = x.*scale;
    y = y.*scale;
    freqs_x = repmat(data.(groups{1}).avg.x_x.freqs,[Nblock,1])';
    freqs_y = repmat(data.(groups{1}).avg.y_y.freqs,[Nblock,1])';
    freqs_x = freqs_x + x;
    freqs_y = freqs_y + y;
    
%% graphs x_x and y_y Bode plots
    for i = 1:length(groups)
        figure('Name',titles{i},'NumberTitle','off')
        subplot(2,2,1) %x input, x output
        errorbar(freqs_x(:,1), 20*log10(data.(groups{i}).avg.x_x.amplitude(1,:)'), data.(groups{i}).avg.x_x.amp_err(1,:,1)',data.(groups{i}).avg.x_x.amp_err(1,:,2)','LineWidth',1.5)
        hold on
        errorbar(freqs_y(:,2:5), 20*log10(data.(groups{i}).avg.y_y.amplitude(2:5,:)'), data.(groups{i}).avg.y_y.amp_err(2:5,:,1)',data.(groups{i}).avg.y_y.amp_err(2:5,:,2)','LineWidth',1.5) 
        errorbar(freqs_x(:,6), 20*log10(data.(groups{i}).avg.x_x.amplitude(6,:)'), data.(groups{i}).avg.x_x.amp_err(6,:,1)',data.(groups{i}).avg.x_x.amp_err(6,:,2)','LineWidth',1.5)
        title('Magnitude (X_{target} -> X_{cursor})'); ylabel('Magnitude (dB)'); xlabel('Frequency (Hz)')
        grid on; set(gca, 'Xscale', 'log')
        xlim([0.09 2.3])
        ylim([-45 0])
      
        subplot(2,2,3)
        errorbar(freqs_x(:,1), data.(groups{i}).avg.x_x.phase(1,:)'*(180/pi), data.(groups{i}).avg.x_x.phase_err(1,:,1)',data.(groups{i}).avg.x_x.phase_err(1,:,2)','LineWidth',1.5)
        hold on
        errorbar(freqs_y(:,2:5), data.(groups{i}).avg.y_y.phase(2:5,:)'*(180/pi), data.(groups{i}).avg.y_y.phase_err(2:5,:,1)',data.(groups{i}).avg.y_y.phase_err(2:5,:,2)','LineWidth',1.5)
        errorbar(freqs_x(:,6), data.(groups{i}).avg.x_x.phase(6,:)'*(180/pi), data.(groups{i}).avg.x_x.phase_err(6,:,1)',data.(groups{i}).avg.x_x.phase_err(5,:,2)','LineWidth',1.5)
        title('Phase (X_{target} -> X_{cursor})'); ylabel('Phase (degrees)'); xlabel('Frequency (Hz)')
        grid on; set(gca, 'Xscale', 'log')
        xlim([0.09 2.3])
        ylim([-450 50])

        subplot(2,2,2) %y input, y output
        errorbar(freqs_y(:,1), 20*log10(data.(groups{i}).avg.y_y.amplitude(1,:)'), data.(groups{i}).avg.y_y.amp_err(1,:,1)',data.(groups{i}).avg.y_y.amp_err(1,:,2)','LineWidth',1.5)
        hold on
        errorbar(freqs_x(:,2:5), 20*log10(data.(groups{i}).avg.x_x.amplitude(2:5,:)'), data.(groups{i}).avg.x_x.amp_err(2:5,:,1)',data.(groups{i}).avg.x_x.amp_err(2:5,:,2)','LineWidth',1.5)
        errorbar(freqs_y(:,6), 20*log10(data.(groups{i}).avg.y_y.amplitude(6,:)'), data.(groups{i}).avg.y_y.amp_err(6,:,1)',data.(groups{i}).avg.y_y.amp_err(6,:,2)','LineWidth',1.5)
        title('Magnitude (Y_{target} -> Y_{cursor})'); ylabel('Magnitude (dB)'); xlabel('Frequency (Hz)')
        grid on; set(gca, 'Xscale', 'log')
        xlim([0.09 2.3])
        ylim([-45 0])

        subplot(2,2,4);
        errorbar(freqs_y(:,1), data.(groups{i}).avg.y_y.phase(1,:)'*(180/pi), data.(groups{i}).avg.y_y.phase_err(1,:,1)',data.(groups{i}).avg.y_y.phase_err(1,:,2)','LineWidth',1.5) 
        hold on
        errorbar(freqs_x(:,2:5), data.(groups{i}).avg.x_x.phase(2:5,:)'*(180/pi), data.(groups{i}).avg.x_x.phase_err(2:5,:,1)',data.(groups{i}).avg.x_x.phase_err(2:5,:,2)','LineWidth',1.5) 
        errorbar(freqs_y(:,6), data.(groups{i}).avg.y_y.phase(6,:)'*(180/pi), data.(groups{i}).avg.y_y.phase_err(6,:,1)',data.(groups{i}).avg.y_y.phase_err(6,:,2)','LineWidth',1.5) 
        title('Phase (Y_{target} -> Y_{cursor})'); ylabel('Phase (degrees)'); xlabel('Frequency (Hz)')
        grid on; set(gca, 'Xscale', 'log')
        xlim([0.09 2.3])
        ylim([-450 50])
        legend(graph_name,'Position',[0 0.4 0.1 0.2])
    end

%% graph aftereffects Bode plots
%     figure('Name','Aftereffects','NumberTitle','off'); 
%     for i = 1:length(groups)
%         subplot(2,2,i); %x input, y output
%         errorbar(freqs_x(:,[1 6]), 20*log10(data.(groups{i}).avg.x_y.amplitude([1 6],:)'), data.(groups{i}).avg.x_y.amp_err([1 6],:,1)',data.(groups{i}).avg.x_y.amp_err([1 6],:,2)','LineWidth',1.5); 
%         title([titles{i},': X_{target} -> Y_{cursor}'],'FontSize',15); ylabel('Magnitude (dB)','FontSize',15); xlabel('Frequency (Hz)','FontSize',15);
%         grid on; set(gca, 'Xscale', 'log');
%         xlim([0.09 2.3]);
%         ylim([-60 0]);
% 
%         subplot(2,2,i+2); %y input, x output
%         errorbar(freqs_y(:,[1 6]), 20*log10(data.(groups{i}).avg.y_x.amplitude([1 6],:)'), data.(groups{i}).avg.y_x.amp_err([1 6],:,1)',data.(groups{i}).avg.y_x.amp_err([1 6],:,2)','LineWidth',1.5); 
%         title([titles{i},': Y_{target} -> X_{cursor}'],'FontSize',15); ylabel('Magnitude (dB)','FontSize',15); xlabel('Frequency (Hz)','FontSize',15);
%         grid on; set(gca, 'Xscale', 'log');
%         xlim([0.09 2.3]);
%         ylim([-60 0]); 
%     end
%     leg = {'Baseline','Aftereffects'};
%     legend(leg,'Position',[0 0.4 0.1 0.1],'FontSize',12);

%% graph hand Bode plots
%     for i = 1:length(groups)
%         figure('Name',titles{i},'NumberTitle','off'); 
%         subplot(2,2,1); hold on;
%         errorbar(freqs_x(:,1), 20*log10(data.(groups{i}).avg.x_x.amplitude(1,:)'), data.(groups{i}).avg.x_x.amp_err(1,:,1)',data.(groups{i}).avg.x_x.amp_err(1,:,2)','LineWidth',1.5);
%         errorbar(freqs_x(:,2:5), 20*log10(data.(groups{i}).avg.x_y.amplitude(2:5,:)'), data.(groups{i}).avg.x_y.amp_err(2:5,:,1)',data.(groups{i}).avg.x_y.amp_err(2:5,:,2)','LineWidth',1.5);
%         errorbar(freqs_x(:,6), 20*log10(data.(groups{i}).avg.x_x.amplitude(6,:)'), data.(groups{i}).avg.x_x.amp_err(6,:,1)',data.(groups{i}).avg.x_x.amp_err(6,:,2)','LineWidth',1.5);
%         title('Magnitude (X_{target} -> X_{hand})','FontSize',15); ylabel('Magnitude (dB)','FontSize',15); xlabel('Frequency (Hz)','FontSize',15);
%         grid on; set(gca, 'Xscale', 'log');
%         xlim([0.09 2.3]);
%         ylim([-40 0]);
%         
%         subplot(2,2,3); hold on;
%         errorbar(freqs_x(:,1), data.(groups{i}).avg.x_x.phase(1,:)'*(180/pi), data.(groups{i}).avg.x_x.phase_err(1,:,1)',data.(groups{i}).avg.x_x.phase_err(1,:,2)','LineWidth',1.5);
%         errorbar(freqs_x(:,2:5), data.(groups{i}).avg.x_y.phase(2:5,:)'*(180/pi), data.(groups{i}).avg.x_y.phase_err(2:5,:,1)',data.(groups{i}).avg.x_y.phase_err(2:5,:,2)','LineWidth',1.5);
%         errorbar(freqs_x(:,6), data.(groups{i}).avg.x_x.phase(6,:)'*(180/pi), data.(groups{i}).avg.x_x.phase_err(6,:,1)',data.(groups{i}).avg.x_x.phase_err(6,:,2)','LineWidth',1.5);
%         title('Phase (X_{target} -> X_{hand})','FontSize',15); ylabel('Phase (degrees)','FontSize',15); xlabel('Frequency (Hz)','FontSize',15);
%         grid on; set(gca, 'Xscale', 'log');
%         xlim([0.09 2.3]);
%         ylim([-800 200]);
%         
%         subplot(2,2,2); hold on;
%         errorbar(freqs_y(:,1), 20*log10(data.(groups{i}).avg.y_y.amplitude(1,:)'), data.(groups{i}).avg.x_x.amp_err(1,:,1)',data.(groups{i}).avg.x_x.amp_err(1,:,2)','LineWidth',1.5);
%         errorbar(freqs_y(:,2:5), 20*log10(data.(groups{i}).avg.y_x.amplitude(2:5,:)'), data.(groups{i}).avg.x_y.amp_err(2:5,:,1)',data.(groups{i}).avg.x_y.amp_err(2:5,:,2)','LineWidth',1.5);
%         errorbar(freqs_y(:,6), 20*log10(data.(groups{i}).avg.y_y.amplitude(6,:)'), data.(groups{i}).avg.x_x.amp_err(6,:,1)',data.(groups{i}).avg.x_x.amp_err(6,:,2)','LineWidth',1.5);
%         title('Magnitude (Y_{target} -> Y_{hand})','FontSize',15); ylabel('Magnitude (dB)','FontSize',15); xlabel('Frequency (Hz)','FontSize',15);
%         grid on; set(gca, 'Xscale', 'log');
%         xlim([0.09 2.3]);
%         ylim([-40 0]);
%         
%         subplot(2,2,4); hold on;
%         errorbar(freqs_y(:,1), data.(groups{i}).avg.y_y.phase(1,:)'*(180/pi), data.(groups{i}).avg.x_x.phase_err(1,:,1)',data.(groups{i}).avg.x_x.phase_err(1,:,2)','LineWidth',1.5);
%         errorbar(freqs_y(:,2:5), data.(groups{i}).avg.y_x.phase(2:5,:)'*(180/pi), data.(groups{i}).avg.x_y.phase_err(2:5,:,1)',data.(groups{i}).avg.x_y.phase_err(2:5,:,2)','LineWidth',1.5);
%         errorbar(freqs_y(:,6), data.(groups{i}).avg.y_y.phase(6,:)'*(180/pi), data.(groups{i}).avg.x_x.phase_err(6,:,1)',data.(groups{i}).avg.x_x.phase_err(6,:,2)','LineWidth',1.5);
%         title('Phase (Y_{target} -> Y_{hand})','FontSize',15); ylabel('Phase (degrees)','FontSize',15); xlabel('Frequency (Hz)','FontSize',15);
%         grid on; set(gca, 'Xscale', 'log');
%         xlim([0.09 2.3]);
%         ylim([-800 200]);
%         legend(graph_name,'Position',[0 0.4 0.1 0.2],'FontSize',12);
%     end
end