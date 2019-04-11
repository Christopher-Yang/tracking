function graph_bode2(data, graph_name, groups)

    rng(1);
    Nblock1 = size(data.(groups{1}).avg.x_x.fft,1);
    Nfreq1 = size(data.(groups{1}).avg.x_x.fft,2);
    a = reshape(datasample([1 -1],Nblock1*Nfreq1),[Nfreq1 Nblock1]);
    x = rand(Nfreq1,Nblock1).*a*0.01;
    y = rand(Nfreq1,Nblock1).*a*0.01;
    scale = repmat(1:Nfreq1,[Nblock1,1])';
    x = x.*scale;
    y = y.*scale;
    freqs_x1 = repmat(data.(groups{1}).avg.x_x.d.freqs,[Nblock1,1])';
    freqs_y1 = repmat(data.(groups{1}).avg.y_y.d.freqs,[Nblock1,1])';
    freqs_x1 = freqs_x1 + x;
    freqs_y1 = freqs_y1 + y;
    
    % for plotting rotated sum of sines in X and Y
    Nblock2 = size(data.(groups{2}).avg.x_x.fft,1);
    Nfreq2 = size(data.(groups{2}).avg.x_x.fft,2);
    a = reshape(datasample([1 -1],Nblock2*Nfreq2),[Nfreq2 Nblock2]);
    x = rand(Nfreq2,Nblock2).*a*0.01;
    y = rand(Nfreq2,Nblock2).*a*0.01;
    scale = repmat(1:Nfreq2,[Nblock2,1])';
    x = x.*scale;
    y = y.*scale;
    freqs_x2 = repmat(data.(groups{2}).avg.x_x.d.freqs,[Nblock2,1])';
    freqs_y2 = repmat(data.(groups{2}).avg.y_y.d.freqs,[Nblock2,1])';
    freqs_x2 = freqs_x2 + x;
    freqs_y2 = freqs_y2 + y;
    
%%    
    figure('Name','Bode plots','NumberTitle','off'); 
    subplot(2,2,1); %x input, x output
    errorbar(freqs_x1, 20*log10(data.(groups{1}).avg.x_x.d.amplitude'), data.(groups{1}).avg.x_x.d.amp_err','LineWidth',1.5); hold on;
%     errorbar(freqs_x2, 20*log10(data.(groups{2}).avg.x_x.d.amplitude'), data.(groups{2}).avg.x_x.d.amp_err','LineWidth',1.5);
    title('Magnitude (X_{target} -> X_{cursor})'); ylabel('Magnitude (dB)'); xlabel('Frequency (Hz)');
    grid on; set(gca, 'Xscale', 'log');
    xlim([0.09 2.3]);
    ylim([-20 5]);
      
    subplot(2,2,3);
    errorbar(freqs_x1, data.(groups{1}).avg.x_x.d.phase'*(180/pi), data.(groups{1}).avg.x_x.d.phase_err','LineWidth',1.5); hold on;
%     errorbar(freqs_x2, data.(groups{2}).avg.x_x.d.phase'*(180/pi), data.(groups{2}).avg.x_x.d.phase_err','LineWidth',1.5);
    title('Phase (X_{target} -> X_{cursor})'); ylabel('Phase (degrees)'); xlabel('Frequency (Hz)');
    grid on; set(gca, 'Xscale', 'log');
    xlim([0.09 2.3]);
    ylim([-450 50]);

    subplot(2,2,2); %y input, y output
    errorbar(freqs_y1, 20*log10(data.(groups{1}).avg.y_y.d.amplitude'), data.(groups{1}).avg.y_y.d.amp_err','LineWidth',1.5); hold on;
%     errorbar(freqs_y2, 20*log10(data.(groups{2}).avg.y_y.d.amplitude'), data.(groups{2}).avg.y_y.d.amp_err','LineWidth',1.5);
    title('Magnitude (Y_{target} -> Y_{cursor})'); ylabel('Magnitude (dB)'); xlabel('Frequency (Hz)');
    grid on; set(gca, 'Xscale', 'log');
    xlim([0.09 2.3]);
    ylim([-20 5]);

    subplot(2,2,4);
    errorbar(freqs_y1, data.(groups{1}).avg.y_y.d.phase'*(180/pi), data.(groups{1}).avg.y_y.d.phase_err','LineWidth',1.5); hold on;
%     errorbar(freqs_y2, data.(groups{2}).avg.y_y.d.phase'*(180/pi), data.(groups{2}).avg.y_y.d.phase_err','LineWidth',1.5);
    title('Phase (Y_{target} -> Y_{cursor})'); ylabel('Phase (degrees)'); xlabel('Frequency (Hz)');
    grid on; set(gca, 'Xscale', 'log');
    xlim([0.09 2.3]);
    ylim([-450 50]);
    legend(graph_name,'Position',[0 0.4 0.1 0.2]);
    
%%
    figure;
    subplot(2,2,1);
    errorbar(freqs_x2, 20*log10(data.(groups{2}).avg.x_x.d.amplitude'), data.(groups{2}).avg.x_x.d.amp_err', 'LineWidth',1.5)
    title('Magnitude (Axis 1)'); ylabel('Magnitude (dB)'); xlabel('Frequency (Hz)');
    grid on; set(gca, 'Xscale', 'log');
    xlim([0.09 2.3]);
    ylim([-20 5]);
    
    subplot(2,2,3);
    errorbar(freqs_x2, data.(groups{2}).avg.x_x.d.phase'*(180/pi), data.(groups{2}).avg.x_x.d.phase_err','LineWidth',1.5);
    title('Phase (Axis 1)'); ylabel('Phase (degrees)'); xlabel('Frequency (Hz)');
    grid on; set(gca, 'Xscale', 'log');
    xlim([0.09 2.3]);
    ylim([-450 50]);
    
    subplot(2,2,2)
    errorbar(freqs_x2, 20*log10(data.(groups{2}).avg.y_y.d.amplitude'), data.(groups{2}).avg.y_y.d.amp_err', 'LineWidth',1.5)
    title('Magnitude (Axis 2)'); ylabel('Magnitude (dB)'); xlabel('Frequency (Hz)');
    grid on; set(gca, 'Xscale', 'log');
    xlim([0.09 2.3]);
    ylim([-20 5]);
    
    subplot(2,2,4);
    errorbar(freqs_x2, data.(groups{2}).avg.y_y.d.phase'*(180/pi), data.(groups{2}).avg.y_y.d.phase_err','LineWidth',1.5);
    title('Phase (Axis 2)'); ylabel('Phase (degrees)'); xlabel('Frequency (Hz)');
    grid on; set(gca, 'Xscale', 'log');
    xlim([0.09 2.3]);
    ylim([-450 50]);
    legend('Rotated 1','Rotated 2');
end