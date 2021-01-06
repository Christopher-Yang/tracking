
% set variables for analysis
rng(3);
groups = {'rot','mir'}; % names of the groups
Ngroup = length(groups); % number of groups
Nsubj = 10; % number of subjects
Nfreq = 7; % number of frequencies
Ntrials = 8; % number of trials
Nsims = 100; % number of simulations to run per observation time
Nreps = 2; % number of base periods to track
delt = 1/130; % simulation time step

% use the same frequencies and amplitudes of sines for simulation
freq = data1.rot{1}.baseline.freqX';
amp = data1.rot{1}.baseline.ampX';

%% compute phase lags from late learning of the main experiment
lag = NaN(Nfreq,Ntrials,Ngroup,Nsubj); % preallocate variable for lag
for j = 1:Ngroup
    for i = 1:Nsubj
        phase = unwrap(data1.(groups{j}){i}.pert4.cursor.phasors.x_x.phase,[],1); % store phases
        lag(:,:,j,i) = -phase./(2*pi*freq) - 0.1; % convert phase from radians to seconds
    end
end
lag = lag - 0.1; % subtract off 0.1 secs off lag to account for system delay

% compute standard deviation, minimum, and maximum lag 
lag_sd = std(lag,[],4);
lag_mu = mean(mean(mean(lag,2),1),4);

%% create target trajectory to track
Ts = 1/(min(freq(:))/2); % base period of the trajectory
Nstep = ceil((Nreps*Ts+6)/delt); % number of simulation time steps
target = NaN(Nsims,Nstep); % preallocate variable for trajectory

% create one target trajectory for each simulation
for i = 1:Nsims
    ph = 2*pi*rand(length(freq),1)-pi; % randomize phases of the sines
    target2 = repmat(amp,[1 Nstep]).*cos(freq*2*pi*(0:delt:Nstep*delt-delt) + repmat(ph,[1 Nstep])); % generate sines
    target(i,:) = squeeze(sum(target2,1)); % sum sines together
end

%% specify infinite-horizon, discrete-time LQR 
% Single joint reaching movements:
G = .14; % dissipative torque
M = 0.1; % moment of inertia
tau = 0.06; % muscle time constant

% create state space model in continuous time
Ac = [0 1 0 0 0
    0 -G/M 1/M 0 0
    0 0 -1/tau 0 0
    0 0 0 0 1
    0 0 0 0 0];
Bc = [0 0 1/tau 0 0]';

mat_c = [Ac Bc;
        0 0 0 0 0 0]; % combine Ac and Bc into one matrix
mat_d = expm(mat_c*delt); % discretize mat_c

% extract discretized A and B matrices for mat_d
A = mat_d(1:5,1:5);
B = mat_d(1:5,end);

% accuracy and effort costs
R = 0.0001; % effort cost
Q = [1 0 0 -1 0
    0 0 0 0 0
    0 0 0 0 0
    -1 0 0 1 0
    0 0 0 0 0]; % accuracy cost

order = size(A,1); % order of the system

% calculate feedback gain matrix, L, by iterating Riccati equation
n = 5000; % number of times to iterate equation
P = zeros(order,order,n);
P(:,:,1) = rand(order);
for i = 2:n % loop for iteraction
    P(:,:,i) =  A'*P(:,:,i-1)*A - (A'*P(:,:,i-1)*B)*inv(R + B'*P(:,:,i-1)*B)*(B'*P(:,:,i-1)*A) + Q;
end
L = inv(R + B'*P(:,:,i)*B)*(B'*P(:,:,i)*A)*.2; % constant feedback gain matrix

%% simulate sum-of-sines tracking
Mu = (0.1:0.01:0.9)'; % intermittent observation times to simulate
Ndelay = length(Mu); % number of observation times

% add normally distributed noise to observation times
noise = normrnd(0,min(min(min(lag_sd))),[Ndelay Nsims]);
Mu2 = repmat(Mu,[1 Nsims]);
Mu2 = Mu2 + noise; % Mu2 is the noised version of Mu
Mu2(Mu2<0) = 0;

% preallocate variables for simulation
coherence = NaN(Ndelay,Nfreq,Nsims); % coherence
x = zeros(order,Nstep); % state vector
u = zeros(Nstep,1); % movement commands
hand = NaN(Nsims,Nstep); % matrix for storing hand position across multiple simulations

disp('Simulating observation time...')
% simulation loop
for k = 1:Ndelay % loop over all observation times
    disp(['   ' num2str(Mu(k)) ' secs'])
    for j = 1:Nsims % repeat simulation Nsims times
        mu = Mu2(k,j); % observation time to be used for current simulation
        intervalTimer = 0; % tracks how much time has elapsed since last observation
        idx = 1; % tracks index of last observation
        
        % initialize hand and target position
        x(1,1) = 0; % hand position
        x(4,1) = target(j,1); % target position
        
        % simulate trajectory for every time step
        for i = 2:Nstep
            u(i) = -L*x(:,i-1); % compute motor command
            x(:,i) = A*x(:,i-1) + B*u(i); % compute next state vector
            
            intervalTimer = intervalTimer + delt; % add elapsed time
            
            % check whether elapsed time has exceeded the observation time
            if intervalTimer > mu % if elapsed time is higher
                intervalTimer = 0; % reset timer
                x(4,i) = target(j,i); % observe a new target position
                idx = i; % record index of observed position
            else
                x(4,i) = target(j,idx); % do not observe a new target position
            end
        end
        
        hand(j,:,k) = x(1,:); % record the hand position for current simulation
        targObserve(j,:,k) = x(4,:);
    end
    
    % compute coherence for each simulation
    for i = 1:Nsims
        coherence(k,:,i) = mscohere(hand(i,6/delt+1:end,k),target(i,6/delt+1:end),blackmanharris(round(Nstep/5)),[],freq,1/delt);
    end
end

%% plot simulation results
% set line colors
col = copper;
col = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);
col2 = [220 20 60
        50 205 50]./255;

coherenceMu = mean(coherence,3);

figure; clf; hold on
plot(1000*[lag_mu(1); lag_mu(1)],[0 1],'LineWidth',2,'Color',col2(1,:)) % lag of rotation group
plot(1000*[lag_mu(2); lag_mu(2)],[0 1],'LineWidth',2,'Color',col2(2,:)) % lag of mirror reversal group
for i = 1:Nfreq
    plot(1000*Mu,coherenceMu(:,i),'-','Color',col(i,:),'LineWidth',2)
end
xticks(100:200:900)
yticks(0:0.25:1)
xlabel('Intermittent Observation Time (ms)')
ylabel('Coherence')
axis([100 900 0 1])
set(gca,'TickDir','out')

%%
f = 81;
trial = 10;

hand_traj = squeeze(hand(:,781:end,:));
target_traj = target(:,781:end)';

hand_fft = fft(hand_traj - repmat(mean(hand_traj,2),[1 5200 1]),[],2);
target_fft = fft(target_traj - repmat(mean(target_traj,1),[5200 1]),[],1);

idx = find(abs(target_fft(:,1))>1.5);
idx = idx(1:end/2);

target_fft = permute(target_fft,[2 1]);
ratio = hand_fft(:,idx,:)./repmat(target_fft(:,idx,:),[1 1 81]);
phase2 = angle(ratio);
lag2 = -phase2./(2*pi*repmat(freq',[Nsims 1 81]));
lag2 = squeeze(mean(mean(lag2,1),2));

figure(20); clf; hold on
plot(Mu*1000,lag2*1000)
ylabel('Lag (ms)')
xlabel('Observation Time (ms)')

%%
n = length(hand_fft);
processed = hand_fft(1:floor(n/2)+1)/n;
processed(2:end-1) = 2*processed(2:end-1);
spectra = abs(processed);

x_axis = 130*(0:n/2)/n;

figure(10); clf; 
subplot(1,2,1); hold on
plot(0:delt:46-delt,hand(1,:,f))
plot(0:delt:46-delt,targObserve(1,:,f),'--k')
plot(0:delt:46-delt,target(1,:),'k')
xlim([10 20])
xlabel('Time (s)')
ylabel('Position (m)')
title(['Observation Time = ' num2str(Mu(f)*1000) ' ms'])
legend({'Hand','Observed target position','True target position'},'Location','southeast')

subplot(1,2,2); hold on
plot(x_axis,spectra)
plot(freq,amp,'.k','MarkerSize',20)
xlim([0 2.2])
title(['Observation Time = ' num2str(Mu(f)*1000) ' ms'])
xlabel('Frequence (Hz)')
ylabel('Amplitude')