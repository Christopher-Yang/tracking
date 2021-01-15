% This function simulates an intermittent catch-up strategy. The optional
% argument "simResults" contains the output of the simulation, which is the
% coherence between the LQR and target movement as well as the hand
% trajectories with a 300 ms sampling period. If simResults is provided, it
% considerably reduces computation time. If simResults isn't provided, then
% it will be computed from scratch. You can save the results of the
% simulation by uncommenting line 142.

function LQR(data,simResults)

% set variables for analysis
rng(3);
Nfreq = 7; % number of frequencies
Nsims = 1000; % number of simulations to run per observation time
Nreps = 2; % number of base periods to track
delt = 1/130; % simulation time step

% use the same frequencies and amplitudes of sines for simulation
freq = data.rot{1}.baseline.freqX';
amp = data.rot{1}.baseline.ampX';

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
Mu = (0.1:0.01:0.5)'; % intermittent observation times to simulate

% only perform simulation if coherence isn't provided 
if nargin == 1
    Ndelay = length(Mu); % number of observation times
    
    % add normally distributed noise to observation times
    noise = normrnd(0,0.02,[Ndelay Nsims]);
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
                        
            hand(j,:) = x(1,:); % record the hand position for current simulation
        end
        
        % store hand trajectory for 300 ms observation time
        if k == 21
            simResults.hand = hand;
        end
        
        % compute coherence for each simulation
        for i = 1:Nsims
            coherence(k,:,i) = mscohere(hand(i,6/delt+1:end),target(i,6/delt+1:end),blackmanharris(round(Nstep/5)),[],freq,1/delt);
        end
    end
    
    simResults.coherence = coherence;
    
    % uncomment the line below to save the coherence results from the
    % simulation
%     save simResults simResults
end

% compute DFTs of the LQR's response
% simNumber = 10;

hand_traj = squeeze(simResults.hand(:,781:end))'; % truncate trajectory to 40 secs
n = size(hand_traj,1);
hand_fft = fft(hand_traj - repmat(mean(hand_traj,1),[n 1])); % take FFT
processed = hand_fft(1:floor(n/2)+1,:)/n;
processed(2:end-1,:) = 2*processed(2:end-1,:);
processed = abs(processed); % amplitude spectrum
spectra = mean(processed,2);
%% plot simulation results
% set line colors
col = copper;
col = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);
col2 = [255 193 7]./255;

% x axis for amplitude spectrum
x_axis = 130*(0:n/2)/n;

% average coherence across simulations
coherenceMu = mean(simResults.coherence,3);

% plot amplitude spectra
figure(14); clf
subplot(1,2,1); hold on
plot(freq,amp*100,'d','Color',col2,'MarkerSize',6,'MarkerFaceColor',col2)
plot(x_axis,spectra*100,'k','LineWidth',0.5)
xlim([0 2.4])
xlabel('Frequency (Hz)')
ylabel('Amplitude (cm)')
xticks(0:2)
yticks(0:2)
set(gca,'TickDir','out')

% plot coherence as a function of observation time
subplot(1,2,2); hold on
plot([300 300],[0 1],'--k','LineWidth',1) 
for i = 1:Nfreq
    plot(1000*Mu,coherenceMu(:,i),'-','Color',col(i,:),'LineWidth',1.5)
end
xticks(100:100:500)
yticks(0:0.25:1)
xlabel('Sampling period (ms)')
ylabel('Coherence')
axis([100 500 0 1])
set(gca,'TickDir','out')

end