function [phasors, raw_fft, processed_fft] = fourier(trajs, trialType)

% compute the raw ffts
cursorHand.xFFT = fft(trajs.cursorHand.x - repmat(mean(trajs.cursorHand.x,1), [size(trajs.cursorHand.x,1) 1])); % hand input to cursor
cursorHand.yFFT = fft(trajs.cursorHand.y - repmat(mean(trajs.cursorHand.y,1), [size(trajs.cursorHand.y,1) 1]));
cursorInput.xFFT = fft(trajs.cursorInput.x - repmat(mean(trajs.cursorInput.x,1), [size(trajs.cursorInput.x,1) 1])); % sinusoidal perturbation of cursor
cursorInput.yFFT = fft(trajs.cursorInput.y - repmat(mean(trajs.cursorInput.y,1), [size(trajs.cursorInput.y,1) 1]));
target.xFFT = fft(trajs.target.x - repmat(mean(trajs.target.x,1), [size(trajs.target.x,1) 1])); % target sines
target.yFFT = fft(trajs.target.y - repmat(mean(trajs.target.y,1), [size(trajs.target.y,1) 1]));

% set variables for analysis
Ntrials = length(trialType); % number of trials

% determine how to analyze raw ffts based on trial type
for i = 1:Ntrials
    
    % analysis for trials with target sines
    in.x = target.xFFT(:,i); % input trajectory
    in.y = target.yFFT(:,i);
    out.x = cursorHand.xFFT(:,i); % output trajectory
    out.y = cursorHand.yFFT(:,i);
    flip = 0;
    
    % compute complex ratios
    phasors.xTarg_x{i} = evaluateFFT(in.x, out.x, flip);
    phasors.xTarg_y{i} = evaluateFFT(in.x, out.y, flip);
    phasors.yTarg_x{i} = evaluateFFT(in.y, out.x, flip);
    phasors.yTarg_y{i} = evaluateFFT(in.y, out.y, flip);
    
    % analysis for trials with cursor sines
    in.x = cursorInput.xFFT(:,i);
    in.y = cursorInput.yFFT(:,i);
    flip = 1;
    
    % compute complex ratios
    phasors.xCurs_x{i} = evaluateFFT(in.x, out.x, flip);
    phasors.xCurs_y{i} = evaluateFFT(in.x, out.y, flip);
    phasors.yCurs_x{i} = evaluateFFT(in.y, out.x, flip);
    phasors.yCurs_y{i} = evaluateFFT(in.y, out.y, flip);

end

% store raw ffts
raw_fft.cursorHand = cursorHand;
raw_fft.cursorInput = cursorInput;
raw_fft.target = target;

% set variables for analysis
n = size(cursorHand.xFFT,1); % number of frequencies
fields = {'xFFT','yFFT'}; % fields in data structure to analyze

% compute amplitude spectra
for i = 1:length(fields)
    cursorHand.(fields{i}) = ampSpectra(cursorHand.(fields{i}),n);
    cursorInput.(fields{i}) = ampSpectra(cursorInput.(fields{i}),n);
    target.(fields{i}) = ampSpectra(target.(fields{i}),n);
end

% store amplitude spectra
processed_fft.cursorHand = cursorHand;
processed_fft.cursorInput = cursorInput;
processed_fft.target = target;
end


% compute complex ratio, gain, and phase
function output = evaluateFFT(in, out, flip)
idx = abs(in) > 10;
idx(length(idx)/2+1:end) = 0;

if flip
    output.ratio = -(out(idx)./in(idx));
else
    output.ratio = out(idx)./in(idx);
end
output.gain = abs(output.ratio);
output.phase = angle(output.ratio);
end

% compute amplitude spectra
function output = ampSpectra(spectrum, n)
output = abs(spectrum(1:floor(n/2)+1,:)/n);
output(2:end-1,:) = 2*output(2:end-1,:);
end