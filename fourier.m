function output = fourier(output_traj,input_traj,Nfreq)

    output = struct;
    switch nargin
        case 1 %calculate DFT
            raw_fft = fft(output_traj-mean(output_traj,1));
            n = size(raw_fft,1);
            output = raw_fft(1:floor(n/2)+1,:)/n;
            output(2:end-1,:) = 2*output(2:end-1,:);
        case 3
            output_fft = NaN(size(output_traj));
            input_fft = NaN(size(input_traj));
            for i = 1:size(output_traj,2)
                output_fft(:,i) = fft(output_traj(:,i)-mean(output_traj(:,i)));
                input_fft(:,i) = fft(input_traj(:,i)-mean(input_traj(:,i)));
            end
            idx = find(abs(input_fft(:,1))>1.5);
            if length(idx) ~= Nfreq*2
                error(['Number of frequencies found (',num2str(length(idx)/2),') does not match the number of input frequencies (',num2str(Nfreq),')'])
            end
            idx = idx(1:Nfreq);
            output.ratio = output_fft(idx,:)./input_fft(idx,:);
            output.gain = abs(output.ratio);
            output.phase = angle(output.ratio);
            output.index = idx;
    end
end