rng(1);

realPart = -2:0.2:2;
imaginaryPart = realPart;
[X,Y] = meshgrid(realPart,imaginaryPart);
comp = X + Y*1j;

Nsamples = 50;
nReal = length(realPart);
nImaginary = length(imaginaryPart);
B_cov = NaN(2,2,nImaginary,nReal);
F_cov = NaN(2,2,nImaginary,nReal,nImaginary*nReal);

B = comp;
F = repmat(reshape(comp,[1 1 numel(comp)]),size(comp));

Hud = -B./(1+B);
Hur = (F+B)./(1+B);

r = normrnd(0, 0.05, [nImaginary nReal Nsamples 2]);
n1 = r(:,:,:,1) + r(:,:,:,2)*1j;

r = normrnd(0, 0.05, [nImaginary nReal nImaginary*nReal Nsamples 2]);
n2 = r(:,:,:,:,1) + r(:,:,:,:,2)*1j;
clear r

Hud_noise = Hud + n1;
Hur_noise = Hur + n2;
clear n1 n2
Hur_noise = permute(Hur_noise,[1 2 4 3]);

B_noise = -Hud_noise ./ (Hud_noise + 1);
F_noise = (Hur_noise + Hud_noise)./(1 + Hud_noise);

e1 = abs(B - B_noise);
e2 = abs(F - permute(F_noise,[1 2 4 3]));

Berror = nanmean(e1,3);
Ferror = nanmean(e2,4);

Ferror = reshape(Ferror,[nImaginary nReal nImaginary nReal]);

% tic
% for i = 1:nImaginary
%     for j = 1:nReal
%         B_cov(:,:,i,j) = cov(real(squeeze(B_noise(i,j,:))),imag(squeeze(B_noise(i,j,:))));
%         for k = 1:nImaginary*nReal
%             F_cov(:,:,i,j,k) = cov(real(squeeze(F_noise(i,j,:,k))),imag(squeeze(F_noise(i,j,:,k))));
%         end
%     end
% end
% toc

F_noise = reshape(F_noise,[nImaginary nReal Nsamples nImaginary nReal]);
% F_cov = reshape(F_cov,[2 2 nImaginary nReal nImaginary nReal]);

F2 = reshape(F,[21 21 21 21]);
Hur2 = reshape(Hur,[21 21 21 21]);
F_noise2 = reshape(F_noise,[21 21 Nsamples 21 21]);
Hur_noise2 = reshape(Hur_noise,[21 21 Nsamples 21 21]);

%%
zMax = max([max(max(Berror)) max(max(max(max(Ferror))))]);

figure(1);
subplot(1,2,1)
surf(X,Y,Berror)
title('B')
xlabel('Real')
ylabel('Imaginary')
zlabel('error')
zlim([0 zMax])

imagIdx = 5;
realIdx = 16;

subplot(1,2,2)
surf(X,Y,squeeze(Ferror(imagIdx,realIdx,:,:)))
title(['F (B = ' num2str(comp(imagIdx,realIdx)) ')'])
xlabel('Real')
ylabel('Imaginary')
zlabel('error')
zlim([0 zMax])

%%
imagIdx = 1;
realIdx = 13;
mu = nanmean(B_noise(imagIdx,realIdx,:));

figure(2); clf
subplot(1,2,1); hold on
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
plot(squeeze(B_noise(imagIdx,realIdx,:)),'.','MarkerSize',15)
plot(B(imagIdx,realIdx),'.','MarkerSize',30)
% error_ellipse(B_cov(:,:,imagIdx,realIdx),[real(mu) imag(mu)],'color','r','conf',.90)
title('B')
axis equal

imagIdx2 = 10;
realIdx2 = 10;
mu = nanmean(F_noise(imagIdx,realIdx,:,imagIdx2,realIdx2));

subplot(1,2,2); hold on
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
plot(squeeze(F_noise(imagIdx,realIdx,:,imagIdx2,realIdx2)),'.','MarkerSize',15)
plot(F2(imagIdx,realIdx,imagIdx2,realIdx2),'.','MarkerSize',30)
% error_ellipse(F_cov(:,:,imagIdx,realIdx,imagIdx2,realIdx2),[real(mu) imag(mu)],'color','r','conf',.90)
title('F')
axis equal

%%
threshold = 0.15;
B_lowError = nansum(nansum(Berror < threshold));
F_lowError = squeeze(nansum(nansum(Ferror < threshold,1),2));

figure(3); clf
surf(X,Y,F_lowError)
xlabel('Real')
ylabel('Imaginary')
zlabel(['error < ' num2str(threshold)])

%%
threshold = 0.15;
comp2 = repmat(comp,[1 1 nImaginary nReal]);

B_good = comp(Berror < threshold);
Hud_good = -B_good ./ (1 + B_good);

for i = 1:nImaginary
    for j = 1:nReal
        test = comp(squeeze(Ferror(i,j,:,:)) < threshold);
        F_good{i}{j} = test;
        Hur_good{i}{j} = (test + comp(i,j)) ./ (1 + comp(i,j));
    end
end

n = length(B_good);
Bidx = [16 7];

col1 = [0 0 0];
col2 = [1 0 0];
col = [linspace(col1(1),col2(1),n)' linspace(col1(2),col2(2),n)' linspace(col1(3),col2(3),n)'];

figure(4); clf
subplot(2,2,1); hold on
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
for i = 1:n
    plot(real(Hud_good(i)),imag(Hud_good(i)),'.','MarkerSize',10,'Color',col(i,:))
end
title('Hud')
axis equal

subplot(2,2,3); hold on
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
for i = 1:n
    plot(real(B_good(i)),imag(B_good(i)),'.','MarkerSize',10,'Color',col(i,:))
end
title('B')
axis equal

a = Hur_good{Bidx(1)}{Bidx(2)};
n = length(a);
col = [linspace(col1(1),col2(1),n)' linspace(col1(2),col2(2),n)' linspace(col1(3),col2(3),n)'];

subplot(2,2,2); hold on
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
for i = 1:n
    plot(real(a(i)),imag(a(i)),'.','MarkerSize',10,'Color',col(i,:))
end
title(['Hur (B = ' num2str(comp(Bidx(1),Bidx(2))) ')'])
axis equal

a = F_good{Bidx(1)}{Bidx(2)};

subplot(2,2,4); hold on
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
for i = 1:n
    plot(real(a(i)),imag(a(i)),'.','MarkerSize',10,'Color',col(i,:))
end
title(['F (B = ' num2str(comp(Bidx(1),Bidx(2))) ')'])
axis equal

%%
figure(5); clf; hold on
plot([-1 1],[0 0],'k')
plot([0 0],[-1 1],'k')
axis equal

idx1 = 5;
idx2 = 2;
idx3 = 11;
idx4 = 6;

plot(real(F2(idx1,idx2,idx3,idx4)),imag(F2(idx1,idx2,idx3,idx4)),'.','MarkerSize',30)
plot(squeeze(F_noise2(idx1,idx2,:,idx3,idx4)),'.','MarkerSize',15)
plot(squeeze(Hur_noise2(idx1,idx2,:,idx3,idx4)),'.','MarkerSize',15)