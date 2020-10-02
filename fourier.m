function [phasors, raw_fft, processed_fft] = fourier(trajs, input, x_axis, trialType)

% compute the raw ffts
cursorHand.xFFT = fft(trajs.cursorHand.x - repmat(mean(trajs.cursorHand.x,1), [size(trajs.cursorHand.x,1) 1])); % hand input to cursor
cursorHand.yFFT = fft(trajs.cursorHand.y - repmat(mean(trajs.cursorHand.y,1), [size(trajs.cursorHand.y,1) 1]));
cursorInput.xFFT = fft(trajs.cursorInput.x - repmat(mean(trajs.cursorInput.x,1), [size(trajs.cursorInput.x,1) 1])); % sinusoidal perturbation of cursor
cursorInput.yFFT = fft(trajs.cursorInput.y - repmat(mean(trajs.cursorInput.y,1), [size(trajs.cursorInput.y,1) 1]));
target.xFFT = fft(trajs.target.x - repmat(mean(trajs.target.x,1), [size(trajs.target.x,1) 1])); % target sines
target.yFFT = fft(trajs.target.y - repmat(mean(trajs.target.y,1), [size(trajs.target.y,1) 1]));

% set variables for analysis
Ntrials = length(trialType); % number of trials
index = findMin(x_axis,input); % indices of the desired frequencies

% determine how to analyze raw ffts based on trial type
for i = 1:Ntrials
    
    % analysis for trials with target sines
    if trialType(i) == 1 || trialType(i) >= 3
        in.x = target.xFFT(:,i); % input trajectory
        in.y = target.yFFT(:,i);
        out.x = cursorHand.xFFT(:,i); % output trajectory
        out.y = cursorHand.yFFT(:,i);
        flip = 0;
        
        % compute complex ratios
        phasors.xTarg_x{i} = evaluateFFT(in.x, out.x, index.tX_freq{i}, flip);
        phasors.xTarg_y{i} = evaluateFFT(in.x, out.y, index.tX_freq{i}, flip);
        phasors.yTarg_x{i} = evaluateFFT(in.y, out.x, index.tY_freq{i}, flip);
        phasors.yTarg_y{i} = evaluateFFT(in.y, out.y, index.tY_freq{i}, flip);

    % if no target sines, set cell array to NaN
    else
        phasors.xTarg_x{i} = NaN;
        phasors.xTarg_y{i} = NaN;
        phasors.yTarg_x{i} = NaN;
        phasors.yTarg_y{i} = NaN;
    end
    
    % analysis for trials with cursor sines
    if trialType(i) >= 2
        in.x = cursorInput.xFFT(:,i);
        in.y = cursorInput.yFFT(:,i);
        out.x = cursorHand.xFFT(:,i);
        out.y = cursorHand.yFFT(:,i);
        flip = 1;
        
        % compute complex ratios
        phasors.xCurs_x{i} = evaluateFFT(in.x, out.x, index.cX_freq{i}, flip);
        phasors.xCurs_y{i} = evaluateFFT(in.x, out.y, index.cX_freq{i}, flip);
        phasors.yCurs_x{i} = evaluateFFT(in.y, out.x, index.cY_freq{i}, flip);
        phasors.yCurs_y{i} = evaluateFFT(in.y, out.y, index.cY_freq{i}, flip);
        
    % if no cursor sines, set cell array to NaN
    else
        phasors.xCurs_x{i} = NaN;
        phasors.xCurs_y{i} = NaN;
        phasors.yCurs_x{i} = NaN;
        phasors.yCurs_y{i} = NaN;
    end
end

% store indices of frequencies
phasors.index = index;

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

% finds the index of the desired frequency 
function output = findMin(x_axis, input)
fields = {'tX_freq','tY_freq','cX_freq','cY_freq'};
Ntrials = length(input.tX_freq);
Nfields = length(fields);
for k = 1:Nfields
    for j = 1:Ntrials
        in = input.(fields{k}){j};
        index = NaN(size(in));
        for i = 1:length(in)
            [m, idx] = min(abs(in(i) - x_axis));
            index(i) = idx;
        end
        output.(fields{k}){j} = index;
    end
end
end

% compute complex ratio, gain, and phase
function output = evaluateFFT(in, out, idx, flip)
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
