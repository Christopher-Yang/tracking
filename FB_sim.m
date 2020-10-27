rng(1);
subj = 2;
for i = 1:6
    phasors(:,i) = data{subj}{i}.phasors.xTarg_x{1}.ratio;
    phasors2(:,i) = data{subj}{i}.phasors.xCurs_x{2}.ratio;
end

Hur = mean(phasors,2);
Hud = mean(phasors2,2);

% sweep over different angles at fixed gain
% angles = flip((0:2*pi/200:2*pi)');
% angles2 = flip((0:2*pi/200:2*pi)');
% N = length(angles);
% amp = linspace(0.8,0.2,N)';
% amp2 = 0.1;
% Hud = amp.*cos(angles) + amp.*sin(angles)*1j;
% Hur = amp.*cos(angles2) + amp.*sin(angles2)*1j;

% sweep over different gains at fixed phase
% amps = (0:0.05:10)';
% angle = 7.9*pi/8;
% Hud = amps.*cos(angle) + amps.*sin(angle)*1j;
% Hur = Hud;

Nfreq = size(Hud,1);
% col1 = [1 0 0];
% col2 = [0.1 0 0];
% col = [linspace(col1(1),col2(1),Nfreq)' linspace(col1(2),col2(2),Nfreq)' linspace(col1(3),col2(3),Nfreq)'];
col = copper;
col = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

Nsamples = 1000;
Hur_boot = NaN(Nfreq,Nsamples);
Hud_boot = NaN(Nfreq,Nsamples);

for i = 1:Nsamples
    r = normrnd(0,0.05,[Nfreq 4]);
    n1 = r(:,1) + r(:,2)*1j;
    n2 = r(:,3) + r(:,4)*1j;
    
    Hur_boot(:,i) = Hur + n1;
    Hud_boot(:,i) = Hud + n2;
end

B = -Hud_boot./(Hud_boot + 1);
F = (Hur_boot + Hud_boot)./(1 + Hud_boot);


figure(1); clf
subplot(2,2,1); hold on
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
for i = 1:Nfreq
    plot(Hud_boot(i,:),'.','MarkerSize',10,'Color',col(i,:))
end
title('Hud')
% axis([-1.5 1.5 -1.5 1.5])
axis equal

subplot(2,2,2); hold on
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
for i = 1:Nfreq
    plot(Hur_boot(i,:),'.','MarkerSize',10,'Color',col(i,:))
end
title('Hur')
% axis([-1.5 1.5 -1.5 1.5])
axis equal

subplot(2,2,3); hold on
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
for i = 1:Nfreq
    plot(B(i,:),'.','MarkerSize',10,'Color',col(i,:))
end
title('B')
% axis([-1.5 1.5 -1.5 1.5])
axis equal

subplot(2,2,4); hold on
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
for i = 1:Nfreq
    plot(F(i,:),'.','MarkerSize',10,'Color',col(i,:))
end
title('F')
% axis([-1.5 1.5 -1.5 1.5])
axis equal
