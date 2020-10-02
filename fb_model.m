s = zpk('s');
G = 1/s;
[err.x, model.x] = predictionError(data, G, 'x');
[err.y, model.y] = predictionError(data, G, 'y');
freq.x = data{1}.sineParams.tX_freq{1};
freq.y = data{1}.sineParams.tY_freq{1};

model.x.B_phase = model.x.B_phase - 360;
model.y.B_phase = model.y.B_phase - 360;
%%
figure(1); clf
subplot(1,2,1)
loglog([1 1],[1 1],'r'); hold on
loglog([1 1],[1 1],'b')
loglog(freq.x,squeeze(err.x.fb_ff(:,1,:)),'Color',[1 0 0 0.4])
loglog(freq.x,squeeze(err.x.fb(:,1,:)),'Color',[0 0 1 0.4])
loglog(freq.x,mean(err.x.fb_ff(:,1,:),3),'Color','r','LineWidth',2)
loglog(freq.x,mean(err.x.fb(:,1,:),3),'Color','b','LineWidth',2)
title('X axis')
xlabel('Frequency (Hz)')
ylabel('Prediction error')
ylim([30 6000])
legend({'FF + FB','FB'},'Location','southeast')
box off

subplot(1,2,2)
loglog(freq.y,squeeze(err.y.fb_ff(:,1,:)),'Color',[1 0 0 0.4]); hold on
loglog(freq.y,squeeze(err.y.fb(:,1,:)),'Color',[0 0 1 0.4])
loglog(freq.y,mean(err.y.fb_ff(:,1,:),3),'Color','r','LineWidth',2)
loglog(freq.y,mean(err.y.fb(:,1,:),3),'Color','b','LineWidth',2)
title('Y axis')
xlabel('Frequency (Hz)')
ylim([30 6000])
box off

%%
Ntrials = length(find(data{1}.trialType==1));
Nsubj = length(data);
Nfreq = length(data{1}.sineParams.tX_freq{1});

fields1 = {'x','y'};
fields2 = {'B_mag','B_phase'};
fields3 = {'F_mag','F_phase'};
% fields4 = {'M_mag_inv','M_phase_inv'};

figure(2); clf
for i = 1:length(fields1)
    for k = 1:length(fields2)
%         B = reshape(permute(model.(fields1{i}).(fields2{k}),[2 3 1]),[Ntrials*Nsubj Nfreq]);
%         F = reshape(permute(model.(fields1{i}).(fields3{k}),[2 3 1]),[Ntrials*Nsubj Nfreq]);
%         M_inv = model.(fields1{i}).(fields4{k});

        B = reshape(permute(model.(fields1{i}).(fields2{k}),[1 3 2]),[Nfreq Nsubj*Ntrials]);
        F = reshape(permute(model.(fields1{i}).(fields3{k}),[1 3 2]),[Nfreq Nsubj*Ntrials]);
        
        plotNum = 2*(k-1)+i;
        
        if k == 1
            offset = 0;
        else
            offset = 360;
        end
        
        subplot(2,2,plotNum); hold on
        plot(freq.(fields1{i}),mean(B,2)-offset,'b','LineWidth',3)
        plot(freq.(fields1{i}),mean(F,2),'r','LineWidth',3)
        plot(freq.(fields1{i}),B-offset,'Color',[0 0 1 0.4])
        plot(freq.(fields1{i}),F,'Color',[1 0 0 0.4])
%         plot(freq.(fields1{i}),M_inv,'r')
        
        switch plotNum
            case 1
                set(gca,'Yscale','log')
                title('X movements')
                ylabel('Gain')
                legend({'feedback','feedforward'},'Location','northwest')
                ylim([0.1 10])
            case 2
                set(gca,'Yscale','log')
                title('Y movements')
                ylim([0.1 10])
            case 3
                ylabel('Phase (degrees)')
                ylim([-420 20])
            case 4
                ylim([-420 20])
        end
        set(gca,'Xscale','log')
    end
end

%%
function [err, model] = predictionError(data, G, axis)
s = data{1}.sineParams;
Nfreq = length(s.tX_freq{1});
Nsubj = length(data);

if strcmp(axis,'x')
    fAxis = {'tX_freq','cX_freq'};
    pAxis = {'xTarg_x','xCurs_x'};
    fftAxis = 'xFFT';
elseif strcmp(axis,'y')
    fAxis = {'tY_freq','cY_freq'};
    pAxis = {'yTarg_y','yCurs_y'};
    fftAxis = 'yFFT';
else
    error('Choose "x" or "y" for axis');
end

freqs = s.(fAxis{1}){1};
rads = 2*pi*freqs;

% G_inv = inv(G);
% [M_mag, M_phase] = bode(G, rads);
% [M_mag_inv, M_phase_inv] = bode(G_inv, rads);
% M_phase = M_phase*pi/180;
% M_phase_inv = M_phase_inv*pi/180;

% M = squeeze(M_mag.*(cos(M_phase) + 1j.*sin(M_phase)));
% M_inv = squeeze(M_mag_inv.*(cos(M_phase_inv) + 1j.*sin(M_phase_inv)));

for i = 1:Nsubj
    a = data{i}.phasors;
    b = data{i}.raw_fft;
    for k = 1:2
        Hur(:,k,i) = a.(pAxis{1}){k}.ratio;
        Hud(:,k,i) = a.(pAxis{2}){k+2}.ratio;
        
        Hur2(:,k,i) = reshape(permute([a.(pAxis{1}){k+4}.ratio a.(pAxis{1}){k+6}.ratio],[2 1]),[Nfreq 1]);
        Hud2(:,k,i) = reshape(permute([a.(pAxis{2}){k+6}.ratio a.(pAxis{2}){k+4}.ratio],[2 1]),[Nfreq 1]);
        
        idx1 = [a.index.(fAxis{1}){k+4}; a.index.(fAxis{1}){k+6}];
        idx2 = [a.index.(fAxis{2}){k+4}; a.index.(fAxis{2}){k+6}];
        
        r = [b.target.(fftAxis)(idx1(1,:),k+4) b.target.(fftAxis)(idx1(2,:),k+6)];
        d = [b.cursorInput.(fftAxis)(idx2(1,:),k+6) b.cursorInput.(fftAxis)(idx2(2,:),k+4)];
        R(:,k,i) = reshape(permute(r,[2 1]),[Nfreq 1]);
        D(:,k,i) = reshape(permute(d,[2 1]),[Nfreq 1]);
    end
end

% B = -Hud./(M.*(1+Hud));
% F = (Hur + M_inv.*Hud)./(1 + Hud);
% Hur_fb = B./(1 + B.*M);

B = -Hud./(1+Hud);
F = (Hur + Hud)./(1 + Hud);
Hur_fb = B./(1 + B);

U.whole_pred = Hur.*R + Hud.*D;
U.fb_pred = Hur_fb.*R + Hud.*D;
U.whole = Hur2.*R + Hud2.*D;

B_mag = abs(B);
B_phase = angle(B);
F_mag = abs(F);
F_phase = angle(F);

err.fb_ff = abs(U.whole - U.whole_pred);
err.fb = abs(U.whole - U.fb_pred);
model.B_mag = B_mag;
model.B_phase = unwrap(B_phase,[],1)*180/pi;
model.F_mag = F_mag;
model.F_phase = unwrap(F_phase,[],1)*180/pi;
end