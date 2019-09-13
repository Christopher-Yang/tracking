%% get trajectories from ifft
% clear all
% load dat
dat = data{subj}.(block);
Nstep = length(data{subj}.(block).time);
freq_axis = 130.004*(0:floor(Nstep/2))/Nstep;
threshold = 0.006;
output = 'Rhand';

subj = 2;
block = 'B26';
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

%% compare target and hand trajectories at different frequencies 
% plot actual trajectories
figure(1); clf;
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
figure(2); clf
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

figure(3); clf
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