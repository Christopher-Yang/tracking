function output = fourier(output_traj,input_traj,Nfreq)

    output = struct;
    switch nargin
        case 1 %calculate DFT
            output.fft = fft(output_traj-mean(output_traj));
            amp = output.fft(1:floor(length(output.fft)/2)+1);
%             pow = (1/(130.004*length(output.fft)))*abs(amp).^2;
            amp = abs(amp/length(output.fft));
            amp(2:end-1) = 2*amp(2:end-1);
%             pow(2:end-1) = 2*pow(2:end-1);
            output.amplitude = amp;
%             output.power = pow;
%             output.amplitude = abs(output.fft/length(output.fft));
%             output.amplitude = output.amplitude(1:floor(length(output.fft)/2)+1);
%             output.amplitude(2:end-1) = 2*output.amplitude(2:end-1);
        case 3
            if isvector(output_traj) == 1 && isvector(input_traj) == 1
                output_fft = fft(output_traj-mean(output_traj));  %subtract the mean to normalize baseline
                input_fft = fft(input_traj-mean(input_traj));
                idx = find(abs(input_fft)>10);
                
                % added for noisy data
                if idx(1) < 5
                    idx = [3 8 14 20 30 38 44];
                else
                    idx = [6 12 18 24 32 42 48];
                end
                
%                 if length(idx) ~= Nfreq*2
%                     error(['Number of frequencies found (',num2str(length(idx)/2),') does not match the number of input frequencies (',num2str(Nfreq),')'])
%                 end
%                 idx = idx(1:Nfreq);
                output.ratio = output_fft(idx)./input_fft(idx);
                output.amplitude = abs(output.ratio);
                output.phase = angle(output.ratio);
                output.index = idx;
            else
                output_fft = NaN(size(output_traj));
                input_fft = NaN(size(input_traj));
                for i = 1:size(output_traj,2)
                    output_fft(:,i) = fft(output_traj(:,i)-mean(output_traj(:,i)));
                    input_fft(:,i) = fft(input_traj(:,i)-mean(input_traj(:,i)));
                end
                idx = find(abs(input_fft(:,1))>10);
                if length(idx) ~= Nfreq*2
                    error(['Number of frequencies found (',num2str(length(idx)/2),') does not match the number of input frequencies (',num2str(Nfreq),')'])
                end
                idx = idx(1:Nfreq);
                output.ratio = output_fft(idx,:)./input_fft(idx,:);
                output.amplitude = abs(output.ratio);
                output.phase = angle(output.ratio);
            end
    end
end