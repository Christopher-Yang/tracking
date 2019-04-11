function graph_bode_solo(data, subj_name, block_name, graph_name)

mag = NaN(7,6);
phase = NaN(7,6);

for i = 1:length(subj_name)
    for j = 1:length(block_name)
        mag(:,j) = data.rot_i.(subj_name{i}).(block_name{j}).y_y.amplitude;
        phase(:,j) = data.rot_i.(subj_name{i}).(block_name{j}).y_y.phase;
    end
end
freqs = data.rot_i.(subj_name{1}).no_rot1.freqs_x;
figure;
semilogx(freqs,20*log10(mag)); grid on;
figure;
semilogx(freqs,my_unwrap(unwrap(phase))*180/pi);
% for i = 1:length(subj_name) 
%     figure('Name',subj_name{i},'NumberTitle','off'); subplot(2,2,1); %x input, x output
%     for j = 1:length(block_name)
%         freqs = data.(subj_name{i}).(block_name{j}).freqs_x;
%         A_out = data.(subj_name{i}).(block_name{j}).cursor.x_fft.amp_interp;
%         A_in = data.(subj_name{i}).(block_name{j}).target.x_fft.amp_interp;
%         if j == 1 || j == length(block_name)
%             semilogx(freqs,20*log10(A_out./A_in),'--');
%         else
%             semilogx(freqs,20*log10(A_out./A_in));
%         end
%         set(gca,'Xtick',freqs);
%         hold on;
%     end
%     title('Magnitude (X -> X)'); ylabel('Magnitude (dB)'); xlabel('Frequency (rad/s)');
%     grid on; legend(graph_name,'Location','southwest'); ylim([-40 5]); xlim([0 2.051]);
%     subplot(2,2,3);
%     for j = 1:length(block_name)
%         freqs = data.(subj_name{i}).(block_name{j}).freqs_x;
%         phi_out = data.(subj_name{i}).(block_name{j}).cursor.x_fft.phase_interp;
%         phi_in = data.(subj_name{i}).(block_name{j}).target.x_fft.phase_interp;
%         if j == 1 || j == length(block_name)
%             semilogx(freqs,my_unwrap(unwrap(phi_out - phi_in))*180/pi,'--');
%         else
%             semilogx(freqs,my_unwrap(unwrap(phi_out - phi_in))*180/pi);
%         end
%         set(gca,'Xtick',freqs);
%         hold on;
%     end
%     title('Phase (X -> X)'); ylabel('Phase (degrees)'); xlabel('Frequency (rad/s)');
%     grid on; legend(graph_name,'Location','southwest'); ylim([-360 0]); xlim([0 2.051]);
% 
%     subplot(2,2,2); %y input, y output
%     for j = 1:length(block_name)
%         freqs = data.(subj_name{i}).(block_name{j}).freqs_y;
%         A_out = data.(subj_name{i}).(block_name{j}).cursor.y_fft.amp_interp;
%         A_in = data.(subj_name{i}).(block_name{j}).target.y_fft.amp_interp;
%         if j == 1 || j == length(block_name)
%             semilogx(freqs,20*log10(A_out./A_in),'--');
%         else
%             semilogx(freqs,20*log10(A_out./A_in));
%         end
%         set(gca,'Xtick',freqs);
%         hold on;
%     end
%     title('Magnitude (Y -> Y)'); ylabel('Magnitude (dB)'); xlabel('Frequency (rad/s)');
%     grid on; legend(graph_name,'Location','southwest'); ylim([-40 5]);
%     subplot(2,2,4);
%     for j = 1:length(block_name)
%         freqs = data.(subj_name{i}).(block_name{j}).freqs_y;
%         phi_out = data.(subj_name{i}).(block_name{j}).cursor.y_fft.phase_interp;
%         phi_in = data.(subj_name{i}).(block_name{j}).target.y_fft.phase_interp;
%         if j == 1 || j == length(block_name)
%             semilogx(freqs,my_unwrap(unwrap(phi_out - phi_in))*180/pi,'--');
%         else
%             semilogx(freqs,my_unwrap(unwrap(phi_out - phi_in))*180/pi);
%         end
%         set(gca,'Xtick',freqs);
%         hold on;
%     end
%     title('Phase (Y -> Y)'); ylabel('Phase (degrees)'); xlabel('Frequency (rad/s)');
%     grid on; legend(graph_name,'Location','southwest'); ylim([-360 0]);
% end

%%

% for i = 1:length(subj_name) 
%     figure('Name',subj_name{i},'NumberTitle','off'); subplot(2,2,1); %x input, y output
%     for j = 1:length(block_name)
%         freqs = data.(subj_name{i}).(block_name{j}).freqs_x;
%         A_out = data.(subj_name{i}).(block_name{j}).cursor.y_fft.amp_interp_opp;
%         A_in = data.(subj_name{i}).(block_name{j}).target.x_fft.amp_interp;
%         if j == 1 || j == length(block_name)
%             semilogx(freqs,20*log10(A_out./A_in),'--');
%         else
%             semilogx(freqs,20*log10(A_out./A_in));
%         end
%         set(gca,'Xtick',freqs);
%         hold on;
%     end
%     title('Magnitude (X -> Y)'); ylabel('Magnitude (dB)'); xlabel('Frequency (rad/s)');
%     grid on; legend(graph_name,'Location','southeast','FontSize',7); ylim([-40 5]);
%     subplot(2,2,3);
%     for j = 1:length(block_name)
%         freqs = data.(subj_name{i}).(block_name{j}).freqs_x;
%         phi_out = data.(subj_name{i}).(block_name{j}).cursor.y_fft.phase_interp_opp;
%         phi_in = data.(subj_name{i}).(block_name{j}).target.x_fft.phase_interp;
%         if j == 1 || j == length(block_name)
%             semilogx(freqs,my_unwrap(unwrap(phi_out - phi_in))*180/pi,'--');
%         else
%             semilogx(freqs,my_unwrap(unwrap(phi_out - phi_in))*180/pi);
%         end
%         set(gca,'Xtick',freqs);
%         hold on;
%     end    
%     title('Phase (X -> Y)'); ylabel('Phase (degrees)'); xlabel('Frequency (rad/s)');
%     grid on; legend(graph_name,'Location','southwest'); ylim([-360 0]);
% 
%     subplot(2,2,2); %y input, x output
%     for j = 1:length(block_name)
%         freqs = data.(subj_name{i}).(block_name{j}).freqs_y;
%         A_out = data.(subj_name{i}).(block_name{j}).cursor.x_fft.amp_interp_opp;
%         A_in = data.(subj_name{i}).(block_name{j}).target.y_fft.amp_interp;
%         if j == 1 || j == length(block_name)
%             semilogx(freqs,20*log10(A_out./A_in),'--');
%         else
%             semilogx(freqs,20*log10(A_out./A_in));
%         end
%         set(gca,'Xtick',freqs); 
%         hold on;
%     end
%     title('Magnitude (Y -> X)'); ylabel('Magnitude (dB)'); xlabel('Frequency (rad/s)');
%     grid on; legend(graph_name,'Location','south','FontSize', 7); ylim([-40 5]);
%     subplot(2,2,4);
%     for j = 1:length(block_name)
%         freqs = data.(subj_name{i}).(block_name{j}).freqs_y;
%         phi_out = data.(subj_name{i}).(block_name{j}).cursor.x_fft.phase_interp_opp;
%         phi_in = data.(subj_name{i}).(block_name{j}).target.y_fft.phase_interp;
%         if j == 1 || j == length(block_name)
%             semilogx(freqs,my_unwrap(unwrap(phi_out - phi_in))*180/pi,'--');
%         else
%             semilogx(freqs,my_unwrap(unwrap(phi_out - phi_in))*180/pi);
%         end
%         set(gca,'Xtick',freqs);
%         hold on;
%     end
%     title('Phase (Y -> X)'); ylabel('Phase (degrees)'); xlabel('Frequency (rad/s)');
%     grid on; legend(graph_name,'Location','southwest'); ylim([-360 0]);
% end
end