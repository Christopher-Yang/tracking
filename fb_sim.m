
interval = pi/16;

angle = -pi:interval:pi-interval;
realPart = cos(angle);
imaginaryPart = sin(angle);
comp1 = permute(realPart + imaginaryPart*1j,[2 1]);
n = length(comp1);
amp = linspace(0,0.7,n)';
comp1 = 0.7.*comp1;

angle = 0:interval:2*pi-interval;
realPart = cos(angle);
imaginaryPart = sin(angle);
comp2 = permute(realPart + imaginaryPart*1j,[2 1]);
n = length(comp2);
amp = linspace(0,0.7,n)';
comp2 = 0.7*comp2;

clear Hud Hur B F
Hud(1,1,:) = comp1;
Hud(2,2,:) = comp1;
Hur(1,1,:) = comp2;
Hur(2,2,:) = comp2;
M = eye(2);

% Hud(1,2,:) = comp1;
% Hud(2,1,:) = comp1;
% Hur(1,2,:) = comp2;
% Hur(2,1,:) = comp2;
% M = [0 1; 1 0];

% Hud(1,2,:) = comp1;
% Hud(2,1,:) = comp1;
% Hur(1,2,:) = comp2;
% Hur(2,1,:) = comp2;
% M = [0 1; -1 0];

% don't fix Hur or Hud
for i = 1:size(Hud,3)
    B(:,:,i) = -Hud(:,:,i)*inv(M*Hud(:,:,i) + eye(2));
    F(:,:,i) = Hur(:,:,i) + B(:,:,i)*(M*Hur(:,:,i) - eye(2));
end


% fix Hur or Hud
% idx = 1;
% fixedHur = false;
% 
% if fixedHur
%     for i = 1:size(Hud,3)
%         B(:,:,i) = -Hud(:,:,i)*inv(M*Hud(:,:,i) + eye(2));
%         F(:,:,i) = Hur(:,:,idx) + B(:,:,i)*(M*Hur(:,:,idx) - eye(2));
%     end
% else
%     for i = 1:size(Hud,3)
%         B(:,:,i) = -Hud(:,:,i)*inv(M*Hud(:,:,i) + eye(2));
%     end
%     for i = 1:size(Hud,3)
%         F(:,:,i) = Hur(:,:,i) + B(:,:,idx)*(M*Hur(:,:,i) - eye(2));
%     end
% end

nReal = numel(realPart);
nImag = numel(imaginaryPart);
col1 = [1 0 0];
col2 = [0 1 0];
col3 = [0 0 1];
col4 = [1 1 0];

map = [linspace(col4(1),col1(1),n)' linspace(col4(2),col1(2),n)' linspace(col4(3),col1(3),n)'];

figure(6); clf
for i = 1:2
    for k = 1:2
        subplot(2,2,2*(i-1)+k); hold on
        plot([-1 1],[0 0],'k')
        plot([0 0],[-1 1],'k')
        for j = 1:n
            plot(real(Hur(i,k,j)),imag(Hur(i,k,j)),'.','MarkerSize',20,'Color',map(j,:))
        end
        title('H_{ur}')
        axis equal
    end
end

figure(7); clf
for i = 1:2
    for k = 1:2
        subplot(2,2,2*(i-1)+k); hold on
        plot([-1 1],[0 0],'k')
        plot([0 0],[-1 1],'k')
        for j = 1:n
            plot(real(Hud(i,k,j)),imag(Hud(i,k,j)),'.','MarkerSize',20,'Color',map(j,:))
        end
        title('H_{ud}')
        axis equal
    end
end

figure(8); clf
for i = 1:2
    for j = 1:2
        subplot(2,2,2*(i-1)+j); hold on
        plot([-1 1],[0 0],'k')
        plot([0 0],[-1 1],'k')
        for k = 1:n
            plot(real(B(i,j,k)),imag(B(i,j,k)),'.','MarkerSize',20,'Color',map(k,:))
        end
        title('B')
        axis equal
    end
end

figure(9); clf
for i = 1:2
    for j = 1:2
        subplot(2,2,2*(i-1)+j); hold on
        plot([-1 1],[0 0],'k')
        plot([0 0],[-1 1],'k')
        for k = 1:n
            plot(real(F(i,j,k)),imag(F(i,j,k)),'.','MarkerSize',20,'Color',map(k,:))
        end
        title('F')
        axis equal
    end
end