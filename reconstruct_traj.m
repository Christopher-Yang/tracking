%% get trajectories from ifft
% clear all
% load dat
subj = 3;
block = 'B26';
dat = data{subj}.(block);
Nstep = length(data{subj}.(block).time);
freq_axis = 130.004*(0:floor(Nstep/2))/Nstep;
threshold = 0.006;
output = 'Rhand';

% construct trajectories on- and off-target frequencies
FT = [dat.target.x_fft.fft dat.target.y_fft.fft];
idx_half1 = abs(FT)>Nstep*threshold/2; % find target indices
idx_half2 = [idx_half1(:,2) idx_half1(:,1)];
idx = sum(idx_half1,2);
idx = logical(repmat(idx,[1 2]));
[t_onFreq, t_offFreq] = rebuild_traj(FT,idx); % on-target frequencies

FT = [dat.(output).x_fft.fft dat.(output).y_fft.fft];
[h_onFreq, h_offFreq] = rebuild_traj(FT,idx); % off-target frequencies
h_half1 = rebuild_traj(FT,idx_half1);
h_half2 = rebuild_traj(FT,idx_half2);

%% get trajectories by pulling amplitude and phase from spectra
tFile = d.(subj).(block).tFile;

% get amplitude spectra
aT = [dat.target.x_fft.amplitude dat.target.y_fft.amplitude];
aH = [dat.(output).x_fft.amplitude dat.(output).y_fft.amplitude];
idxT = aT>threshold;
idxH = sum(idxT,2);
idxH = logical(repmat(idxH,[1 2]));

% calculate phase spectra
FT = [dat.target.x_fft.fft dat.target.y_fft.fft]/Nstep;
FT = FT(1:floor(Nstep/2)+1,:);
FT(~idxT) = 0;
pT = angle(FT);
FT = [dat.(output).x_fft.fft dat.(output).y_fft.fft]/Nstep;
FT = FT(1:floor(Nstep/2)+1,:);
FT(~idxH) = 0;
pH = angle(FT);

% plot amplitude and phase spectra
figure(1); clf
for i = 1:2
    subplot(2,2,i); hold on
    stem(freq_axis,aT(:,i))
    stem(freq_axis,aH(:,i))
    ylabel('Amplitude (m)')
    xlim([0 2.2])

    subplot(2,2,i+2); hold on
    stem(freq_axis,pT(:,i)*180/pi)
    stem(freq_axis,pH(:,i)*180/pi)
    xlabel('Frequency (Hz)')
    ylabel('Phase (degrees)')
    xlim([0 2.2])
end

% simulate target and hand trajectories
simTime = dat.time'/1000;
freqsT = tFile(1:14)';
freqsH = sort(freqsT);
freqsH = [freqsH; freqsH];

ampT = repmat(aT(idxT),[1 Nstep]);
ampH = repmat(aH(idxH),[1 Nstep]);
phaseT = repmat(pT(idxT),[1 Nstep]);
phaseH = repmat(pH(idxH),[1 Nstep]);

tSines = ampT.*cos(2*pi*freqsT*simTime + phaseT);
target(:,1) = sum(tSines(1:length(freqsT)/2,:));
target(:,2) = sum(tSines(length(freqsT)/2+1:end,:));
hSines = ampH.*cos(2*pi*freqsH*simTime + phaseH);
hand(:,1) = sum(hSines(1:length(freqsH)/2,:));
hand(:,2) = sum(hSines(length(freqsH)/2+1:end,:));

%% compare target and hand trajectories at different frequencies 
% plot actual trajectories
figure(2); clf;
subplot(2,3,1); hold on
plot(dat.target.x_pos,dat.target.y_pos) % actual target trajectory
plot(dat.(output).x_pos,dat.(output).y_pos) % actual output trajectory
xlabel('X position (m)')
ylabel('Y position (m)')
axis square
title('Actual trajectories')
axis([-.1 .1 -.1 .1]) 

% plot on-target frequency trajectories
subplot(2,3,4); hold on
plot(t_onFreq(:,1),t_onFreq(:,2)) % ifft target and hand trajectories
plot(h_onFreq(:,1),h_onFreq(:,2))
% plot(target(:,1),target(:,2)) % simulated target and hand trajectories
% plot(hand(:,1),hand(:,2))
xlabel('X position (m)')
ylabel('Y position (m)')
axis square
title('Movement at all target frequencies')
axis([-.1 .1 -.1 .1]) 

subplot(2,3,2); hold on
plot(t_onFreq(:,1),t_onFreq(:,2)) % ifft target and hand trajectories
plot(h_half1(:,1),h_half1(:,2))
xlabel('X position (m)')
ylabel('Y position (m)')
axis square
title('Movement at baseline frequencies')
axis([-.1 .1 -.1 .1])

subplot(2,3,5); hold on
plot(t_onFreq(:,1),t_onFreq(:,2)) % ifft target and hand trajectories
plot(h_half2(:,1),h_half2(:,2))
xlabel('X position (m)')
ylabel('Y position (m)')
axis square
title('Movement at compensated frequencies')
axis([-.1 .1 -.1 .1])

% plot off-target frequency trajectories
subplot(2,3,[3 6]); hold on
plot(t_onFreq(1:end-50,1)-h_onFreq(51:end,1),t_onFreq(1:end-50,2)-h_onFreq(51:end,2)) % delayed difference between target and hand trajectories
plot(h_offFreq(:,1),h_offFreq(:,2)) % ifft target and hand trajectories
% plot(dat.(output).x_pos - hand(:,1),dat.(output).y_pos - hand(:,2)) % subtract time domain trajectories
xlabel('X position (m)')
ylabel('Y position (m)')
axis square
title('Movement at non-target frequencies')
axis([-.1 .1 -.1 .1]) 

%% plot transformation matrices
figure(3); clf
col1 = [1 0 0];
col2 = [1 1 1];
Nstep = 100;
map1 = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2),Nstep)', linspace(col1(3),col2(3),Nstep)'];

col1 = [1 1 1];
col2 = [0 0 1];
map2 = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2),Nstep)', linspace(col1(3),col2(3),Nstep)'];

map = [map1; map2];
clims = [-1 1];
for i = 1:3
    if i == 1
        in = [dat.target.x_pos dat.target.y_pos]';
        out = [dat.(output).x_pos dat.(output).y_pos]';
    elseif i == 2
        in = t_onFreq';
        out = h_onFreq';
    else
        in = [t_onFreq(1:end-50,1)-h_onFreq(51:end,1) t_onFreq(1:end-50,2)-h_onFreq(51:end,2)]';
        out = [h_offFreq(51:end,1) h_offFreq(51:end,2)]';
    end
    
    paramsInit = [1 0 0 1];
    err = @(params) sim_error(params,out,in);
    [params_opt,fval] = fmincon(err,paramsInit);
    rotMat = [params_opt(1:2); params_opt(3:4)];

    subplot(1,3,i)
    imagesc(rotMat,clims)
    colormap(map)
    axis square
end

%% plot single sinusoids
simTime2 = simTime(1:1000);

figure(4); clf
for i = 1:14
    if i == 1
        subplot(1,2,1); hold on
        m = 1;
    elseif i == 8
        subplot(1,2,2); hold on
        m = 1;
    end
    target2 = 0.5*cos(2*pi*freqsT(i)*simTime2 + phaseT(i));
    hand2 = 0.5*cos(2*pi*freqsH(i)*simTime2 + phaseH(i));
    plot(simTime2,target2+m*1.5,'k')
    plot(simTime2,hand2+m*1.5,'r')
    m = m + 1;
end

%% simulate tracking experiment
simTime = 0:(1/130.004):46; % simulate an extra six seconds of data to match experiment's phase
freqs = tFile(1:14)';
amp = repmat(tFile(15:28)',[1 length(simTime)]);
phase = repmat(tFile(29:42)',[1 length(simTime)]);
Nfreq = length(freqs);

% sines = amp.*cos(2*pi*freqs*simTime+phase);
% target(:,1) = sum(sines(1:7,end-Nstep+1:end));
% target(:,2) = sum(sines(8:14,end-Nstep+1:end));

sines = amp([1 8],:).*cos(2*pi*freqs([1 8])*simTime+phase([1 8],:));
target(:,1) = sines(1,end-Nstep+1:end);
target(:,2) = sines(2,end-Nstep+1:end);

t_fft = fft(target - repmat(mean(target,1),[Nstep 1]))/Nstep;
aT = abs(t_fft(1:floor(length(t_fft)/2)+1,:));
aT(2:end-1,:) = 2*aT(2:end-1,:);
idx = aT>threshold;

t_fft2 = t_fft(1:floor(Nstep/2)+1,:);
t_fft2(~idx) = 0;
pT = angle(t_fft2);

figure(1); clf
subplot(2,1,1)
stem(freq_axis,aT)
xlim([0 2.2])

subplot(2,1,2)
stem(freq_axis,pT)
xlim([0 2.2])

simTime = 0:(1/130.004):40;
ampT = repmat(aT(idx),[1 Nstep]);
phaseT = repmat(pT(idx),[1 Nstep]);

tSines = ampT.*cos(2*pi*freqs([1 8])*simTime + phaseT);
target2(:,1) = tSines(1,:);
target2(:,2) = tSines(2,:);

figure(2); clf; hold on
plot(target(:,1),target(:,2))
plot(target2(:,1),target2(:,2))

function e = sim_error(params,hand,target)
    rotMat = [params(1:2); params(3:4)];
    rotTarget = rotMat*target;
    d = (hand(:,51:end)-rotTarget(:,1:end-50)).^2;
    e = mean(sum(d,1));
end

function [onFreq, offFreq] = rebuild_traj(FT,idx)
    FT_onFreq = FT;
    FT_onFreq(~idx) = 0;
    onFreq = ifft(FT_onFreq);
    FT_offFreq = FT;
    FT_offFreq(idx) = 0;
    offFreq = ifft(FT_offFreq);
end