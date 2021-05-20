rng(3);
Nfreq = 7; % number of frequencies
Nsims = 20; % number of simulations to run per observation time
Nreps = 2; % number of base periods to track
delt = 1/130; % simulation time step

% use the same frequencies and amplitudes of sines for simulation
freq = data1.rot{1}.baseline.freqX';
amp = data1.rot{1}.baseline.ampX';

Ts = 1/(min(freq(:))/2); % base period of the trajectory
Nstep = ceil((Nreps*Ts+6)/delt); % number of simulation time steps
target = NaN(Nsims,Nstep); % preallocate variable for trajectory

% create one target trajectory for each simulation
for i = 1:Nsims
    ph = 2*pi*rand(length(freq),1)-pi; % randomize phases of the sines
    target2 = repmat(amp,[1 Nstep]).*cos(freq*2*pi*(0:delt:Nstep*delt-delt) + repmat(ph,[1 Nstep])); % generate sines
    target(i,:) = squeeze(sum(target2,1)); % sum sines together
end

clear A
delayAll = 1:13;
for j = 1:length(delayAll)
    delay = delayAll(j);
    for i = 1:Nsims
        t = target(i,:);
        pos = t(2:end);
        vel = diff(t)./delt;
        
        state(:,:,i) = [pos; vel];
        state_t1 = state(:,1:end-delay,i);
        state_t2 = state(:,1+delay:end,i);
        
        A(:,:,j,i) = state_t2*pinv(state_t1);
    end
end

A_mu = mean(A,4);
%%
figure(1); clf;
subplot(2,2,1)
plot(delt:delt:delayAll(end)*delt,squeeze(A_mu(1,1,:)))
xlabel('\Delta (s)')
ylim([-3 1.2])
title('pos_{t+\Delta} = a_{11}*pos_t')

subplot(2,2,2)
plot(delt:delt:delayAll(end)*delt,squeeze(A_mu(1,2,:)))
xlabel('\Delta (s)')
ylim([-3 1.2])
title('pos_{t+\Delta} = a_{12}*vel_t')

subplot(2,2,3)
plot(delt:delt:delayAll(end)*delt,squeeze(A_mu(2,1,:)))
xlabel('\Delta (s)')
ylim([-3 1.2])
title('vel_{t+\Delta} = a_{21}*pos_t')

subplot(2,2,4)
plot(delt:delt:delayAll(end)*delt,squeeze(A_mu(2,2,:)))
xlabel('\Delta (s)')
ylim([-3 1.2])
title('vel_{t+\Delta} = a_{22}*vel_t')
%%
Nsteps = size(target,2);

sim = 4;

state_pred = NaN(2,Nsteps);
state_pred(:,1) = state(:,1,sim);

for i = 2:Nsteps
    state_pred(:,i) = A_mu(:,:,end)*state_pred(:,i-1);
end

figure(1); clf; hold on
plot(state(1,:,sim))
plot(state_pred(1,:))
xlim([0 100])