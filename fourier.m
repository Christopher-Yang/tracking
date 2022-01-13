function processed = fourier(output_traj,input_traj,Nfreq)

    output_fft = NaN(size(output_traj));
    input_fft = NaN(size(input_traj));
    for i = 1:size(output_traj,2)
        output_fft(:,i) = fft(output_traj(:,i));
        input_fft(:,i) = fft(input_traj(:,i));
    end
    idx = find(abs(input_fft(:,1))>1.5);
    if length(idx) ~= Nfreq*2
        error(['Number of frequencies found (',num2str(length(idx)/2),') does not match the number of input frequencies (',num2str(Nfreq),')'])
    end
    idx = idx(1:Nfreq);
    processed.ratio = output_fft(idx,:)./input_fft(idx,:);
    processed.gain = abs(processed.ratio);
    processed.phase = unwrap(angle(processed.ratio));
    processed.index = idx;
    
end