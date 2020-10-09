rng(1);

% Hur = data{1}{6}.phasors.xTarg_x{1}.ratio;
% Hud = data{1}{6}.phasors.xCurs_x{2}.ratio;

% sweep over different angles at fixed gain
angles = flip((0:2*pi/200:2*pi)');
angles2 = flip((0:2*pi/200:2*pi)');
N = length(angles);
amp = linspace(0.8,0.2,N)';
amp2 = 0.1;
Hud = amp.*cos(angles) + amp.*sin(angles)*1j;
Hur = amp.*cos(angles2) + amp.*sin(angles2)*1j;

% sweep over different gains at fixed phase
% amps = (0:0.05:10)';
% angle = 7.9*pi/8;
% Hud = amps.*cos(angle) + amps.*sin(angle)*1j;
% Hur = Hud;

Nfreq = size(Hud,1);
col1 = [1 0 0];
col2 = [0.1 0 0];
col = [linspace(col1(1),col2(1),Nfreq)' linspace(col1(2),col2(2),Nfreq)' linspace(col1(3),col2(3),Nfreq)'];

Nsamples = 1000;
Hur_boot = NaN(Nfreq,Nsamples);
Hud_boot = NaN(Nfreq,Nsamples);

for i = 1:Nsamples
    r = normrnd(0,0.01,[Nfreq 4]);
    n1 = r(:,1) + r(:,2)*1j;
    n2 = r(:,3) + r(:,4)*1j;
    
    Hur_boot(:,i) = Hur + n1;
    Hud_boot(:,i) = Hud + n2;
end

B = -Hud_boot./(Hud_boot + 1);
F = (Hur_boot + Hud_boot)./(1 + Hud_boot);


figure(1); clf
subplot(1,2,1); hold on
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
for i = 1:Nfreq
    plot(Hur_boot(i,:),'.','MarkerSize',10,'Color',col(i,:))
end
title('Hur')
axis([-1.5 1.5 -1.5 1.5])
axis square

subplot(1,2,2); hold on
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
for i = 1:Nfreq
    plot(Hud_boot(i,:),'.','MarkerSize',10,'Color',col(i,:))
end
title('Hud')
axis([-1.5 1.5 -1.5 1.5])
axis square

figure(2); clf
subplot(1,2,1); hold on
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
for i = 1:Nfreq
    plot(F(i,:),'.','MarkerSize',10,'Color',col(i,:))
end
title('F')
% axis([-1.5 1.5 -1.5 1.5])
axis equal

subplot(1,2,2); hold on
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
for i = 1:Nfreq
    plot(B(i,:),'.','MarkerSize',10,'Color',col(i,:))
end
title('B')
% axis([-1.5 1.5 -1.5 1.5])
axis equal
