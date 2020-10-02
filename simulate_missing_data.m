
fields2 = {'tX_freq','tY_freq','tX_amp','tY_amp','tX_phase','tY_phase','cX_freq','cY_freq','cX_amp','cY_amp','cX_phase','cY_phase'};

for j = 1:length(fields2)
    input.(fields2{j}) = d{1}.(fields2{j});
end

traj = reconstruct_cursor(input, d{1}.time, d{1}.frameRate);
Nsamples = length(d{1}.time);
time = [.25 .4 .6 .8 1];

for i = 1:length(time)
    start = 1500;
    nanSamples = start:start+time(i)*60;
    traj.y(nanSamples, 3) = NaN;
    traj.y(nanSamples(1)-1:nanSamples(end)+1,3) = linspace(traj.y(nanSamples(1)-1,3),traj.y(nanSamples(end)+1,3),length(nanSamples)+2);

    traj.xFFT = fft(traj.x - repmat(mean(traj.x,1),[Nsamples 1]));
    traj.yFFT = fft(traj.y - repmat(mean(traj.y,1),[Nsamples 1]));

    traj.xAmp = abs(traj.xFFT(1:floor(Nsamples/2)+1,:)/Nsamples);
    traj.xAmp(2:end-1,:) = 2*traj.xAmp(2:end-1,:);
    traj.yAmp = abs(traj.yFFT(1:floor(Nsamples/2)+1,:)/Nsamples);
    traj.yAmp(2:end-1,:) = 2*traj.yAmp(2:end-1,:);
    
    amp(:,i) = traj.yAmp(:,3);
end

figure(1); clf
for i = 1:length(time)
    subplot(1,length(time),i); hold on
    plot(data{1}.x_axis,amp(:,i))
    plot(data{1}.sineParams.cY_freq{3},data{1}.sineParams.cY_amp{3},'o')
    axis([0 2 0 1.65])
    xlabel('Frequency (Hz)')
    title(['Frame drop duration = ' num2str(time(i)*1000) ' ms'])
    if i == 1
        ylabel('Amplitude (cm)')
    end
end

function output = reconstruct_cursor(input, time, frameRate)
    Ntrials = size(time,2);
    outputX = NaN(size(time));
    outputY = NaN(size(time));
    
    for i = 1:Ntrials
        t = 0:1/frameRate:40-(1/frameRate);
        x = repmat(input.cX_amp{i}',[1 length(t)]).*cos(2*pi*input.cX_freq{i}'*t + repmat(input.cX_phase{i}',[1 length(t)]));
        y = repmat(input.cY_amp{i}',[1 length(t)]).*cos(2*pi*input.cY_freq{i}'*t + repmat(input.cY_phase{i}',[1 length(t)]));
        x = sum(x,1);
        y = sum(y,1);

%         idx = find(isnan(time(:,i)));
%         x(idx) = NaN;
%         y(idx) = NaN;
        
        output.x(:,i) = x';
        output.y(:,i) = y';
    end
end