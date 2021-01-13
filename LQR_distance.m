function LQR_distance(data)

% set variables for analysis
rng(10);
Nfreq = 7; % number of frequencies
Nsims = 500; % number of simulations to run per observation time
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
Mu = (0:0.002:0.05)';
Ndelay = length(Mu);

% add normally distributed noise to observation times
noise = normrnd(0,0.002,[Ndelay Nsims]);
Mu2 = repmat(Mu,[1 Nsims]);
Mu2 = Mu2 + noise; % Mu2 is the noised version of Mu
Mu2(Mu2<0) = 0;

% preallocate variables for simulation
coherence = NaN(Ndelay,Nfreq,Nsims); % coherence
pauseObserve = NaN(Nstep,Ndelay,Nsims);
x = zeros(order,Nstep); % state vector
u = zeros(Nstep,1); % movement commands
hand = NaN(Nstep,Ndelay,Nsims); % matrix for storing hand position across multiple simulations
disp('Simulating observation distance...')

for k = 1:Ndelay
    disp(['   ' num2str(Mu(k)*100) ' cm'])
    for j = 1:Nsims
        mu = Mu2(k,j);
        idx = 1;
        
        x(1,1) = 0; % hand position
        x(4,1) = target(j,1); % target position
        pauseObserve(1,k,j) = 0;
        
        % simulate trajectory
        for i = 2:Nstep
            u(i) = -L*x(:,i-1);
            x(:,i) = A*x(:,i-1) + B*u(i);
            
            error = abs(x(1,i) - target(j,i));
            
            if error > mu
                x(4,i) = target(j,i);
                idx = i;
                pauseObserve(i,k,j) = 0;
            else
                x(4,i) = target(j,idx);
                pauseObserve(i,k,j) = 1;
            end
        end
        
        hand(:,k,j) = x(1,:);
    end
    
    for i = 1:Nsims
        coherence(k,:,i) = mscohere(hand(6/delt+1:end,k,i),target(i,6/delt+1:end),blackmanharris(round(Nstep/5)),[],freq,1/delt);
    end
end

%% create look-up table for determining minimum observation distance
for k = 1:Ndelay
    idx = 1;
    for j = 1:Nsims
        timeCount = 0;
        for i = 1:Nstep
            if pauseObserve(i,k,j) == 1 
                timeCount = timeCount + delt;
                if i == Nstep
                    observationTimes{k}(idx) = timeCount;
                    idx = idx + 1;
                end
            else
                if i ~= 1 && pauseObserve(i-1,k,j) == 1
                    observationTimes{k}(idx) = timeCount;
                    timeCount = 0;
                    idx = idx + 1;
                end
            end
        end
    end
    timesMu(k) = mean(observationTimes{k});
end

% find the observation distance that corresponds closest to a 300 ms
% observation time and store it as "idx"
[~, idx] = min(abs(0.3 - timesMu));
mu = Mu(idx);

%% plot results of simulation
% set line colors
col = copper;
col = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

col2 = [50 205 50]./255;

figure(2); clf; hold on
plot([mu mu]*100,[0 1],'--k','LineWidth',1)
for i = 1:Nfreq
    plot(Mu*100,mean(coherence(:,i,:),3),'Color',col(i,:),'LineWidth',1.5)
end
xlabel('Distance threshold (cm)')
ylabel('Coherence')
legend({'Predicted minimum threshold'},'Location','southwest')
yticks(0:0.25:1)
xticks(1:6)

end