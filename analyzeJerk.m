subj = 4;
trial = 6;
block = 'baseline';

delt = 1/130;
t = (0:delt:40-delt)';
[b, a] = butter(3,14/(130/2));

noise = dsp.ColoredNoise('Color','pink','SamplesPerFrame',5200*8);

z = data1.mir{subj}.(block).cursor;
% z = data1.rot{subj}.(block).target;
% x = sgolayfilt(z.x_pos,3,9);
% y = sgolayfilt(z.y_pos,3,9);

x = filter(b, a, z.x_pos);
y = filter(b, a, z.y_pos);

% rng(5);
% x = filter(b, a, rand(5200,8));
% y = filter(b, a, rand(5200,8));
% x = filter(b, a, reshape(noise(),[5200 8]));
% y = filter(b, a, reshape(noise(),[5200 8]));

vel_x = diff(x)./delt;
vel_y = diff(y)./delt;

vel = filter(b, a, sqrt(vel_x.^2 + vel_y.^2));
acc = filter(b, a, diff(vel)./delt);
jerk = filter(b, a, diff(acc)./delt);

for i = 1:8
    [pks, locs] = findpeaks(jerk(:,i));
    idx = pks > 30;
%     peakTime{i} = locs(idx);
    peaks{i} = pks;
    peakTime{i} = locs;
    duration{i} = diff(peakTime{i}*delt);
end

jerk_fft = fft(jerk-repmat(mean(jerk,1),[5197 1]),[],1);
n = size(jerk_fft,1);
processed = jerk_fft(1:floor(n/2)+1,:)/n;
processed(2:end-1,:) = 2*processed(2:end-1,:);

x_axis = 130*(0:n/2)/n;

figure(3); clf
subplot(1,2,1)
histogram(duration{trial},'Normalization','probability')
xlabel('Time b/w peaks')
ylim([0 0.5])

subplot(1,2,2)
histogram(peaks{trial},0:25:1000,'Normalization','probability')
xlabel('Peak amplitude')
ylim([0 0.7])

% figure(1); clf
% subplot(2,4,1); hold on
% plot(x(:,trial),y(:,trial))
% axis equal
% title('Position')
% 
% subplot(2,4,2); hold on
% plot(1:delt:40-2*delt,vel(131:end,trial))
% title('Velocity')
% 
% subplot(2,4,3); hold on
% plot(1:delt:40-3*delt,acc(131:end,trial))
% title('Acceleration')
% 
% subplot(2,4,4); hold on
% plot(1:delt:40-4*delt,jerk(131:end,trial))
% title('Jerk')
% 
% subplot(2,4,5:8)
% plot(x_axis, abs(processed(:,trial)))
% title('FFT')
% xlabel('Frequency (Hz)')

times = t(peakTime{trial})';

figure(2); clf; hold on
subplot(4,1,1); hold on
plot(t(131:end),x(131:end,trial))
title('Position')
xlim([1 21])

subplot(4,1,2); hold on
vel2 = vel(131:end,trial);
% plot([times; times],[0 max(vel2)],'k')
plot(t(131:end-1),vel2,'LineWidth',1.5)
title('Velocity')
xlim([1 21])

subplot(4,1,3); hold on
acc2 = acc(131:end,trial);
% plot([times; times],[min(acc2) max(acc2)],'k')
plot(t(131:end-2),acc2,'LineWidth',1.5)
title('Acceleration')
xlim([1 21])

subplot(4,1,4); hold on
jerk2 = jerk(131:end,trial);
plot([times; times],[min(jerk2) max(jerk2)],'k')
plot(t(131:end-3),jerk2,'LineWidth',1.5)
title('Jerk')
xlim([1 21])
