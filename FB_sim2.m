% realPart = -1.3:0.1:1.3;
% imaginaryPart = realPart;
% % [X,Y] = meshgrid(realPart,imaginaryPart);
% % comp = X + Y*1j;
% comp = permute(realPart + imaginaryPart*1j,[2 1]);
% n = length(comp);

interval = pi/16;
angle = 0:interval:2*pi-interval;
realPart = cos(angle);
imaginaryPart = sin(angle);
comp = permute(realPart + imaginaryPart*1j,[2 1]);
n = length(comp);
amp = linspace(0.5,0.8,n)';
comp = 0.7.*comp;

clear Hud Hur B F
% Hud(1,1,:) = comp;
% Hud(2,2,:) = comp;
% Hur(1,1,:) = comp;
% Hur(2,2,:) = comp;
% M = eye(2);

Hud(1,2,:) = comp;
Hud(2,1,:) = comp;
Hur(1,2,:) = comp;
Hur(2,1,:) = comp;
M = [0 1; 1 0];

idx = 1;
fixedHur = true;

% disturbance after M
if fixedHur
    for i = 1:size(Hud,3)
        B(:,:,i) = -Hud(:,:,i)*inv(M*Hud(:,:,i) + eye(2));
        F(:,:,i) = Hur(:,:,idx) + B(:,:,i)*(M*Hur(:,:,idx) - eye(2));
    end
else
    for i = 1:size(Hud,3)
        B(:,:,i) = -Hud(:,:,i)*inv(M*Hud(:,:,i) + eye(2));
    end
    for i = 1:size(Hud,3)
        F(:,:,i) = Hur(:,:,i) + B(:,:,idx)*(M*Hur(:,:,i) - eye(2));
    end
end

% disturbance before M
% if fixedHur
%     for i = 1:size(Hud,3)
%         B(:,:,i) = -Hud(:,:,i)*inv(M*Hud(:,:,i) + M);
% %         B(:,:,i) = -Hud(:,:,i)*M*inv(Hud(:,:,i) + eye(2));
%         F(:,:,i) = Hur(:,:,idx) + B(:,:,i)*(M*Hur(:,:,idx) - eye(2));
%     end
% else
%     for i = 1:size(Hud,3)
%         B(:,:,i) = -Hud(:,:,i)*inv(M*Hud(:,:,i) + M);
% %         B(:,:,i) = -Hud(:,:,i)*M*inv(Hud(:,:,i) + eye(2));
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

map = [linspace(col1(1),col4(1),n)' linspace(col1(2),col4(2),n)' linspace(col1(3),col4(3),n)'];

ax1 = 1;
ax2 = 1;

figure(7); clf
for i = 1:2
    for k = 1:2
        subplot(2,2,2*(i-1)+k); hold on
        plot([-1 1],[0 0],'k')
        plot([0 0],[-1 1],'k')
        % for i = 1:n
            for j = 1:n
        %         plot(real(comp(i,j)),imag(comp(i,j)),'.','MarkerSize',20,'Color',squeeze(map(i,j,:)))
                plot(real(Hud(i,k,j)),imag(Hud(i,k,j)),'.','MarkerSize',20,'Color',map(j,:))
            end
        % end
        title('H_{ur} and H_{ud}')
        axis equal
    end
end

figure(8); clf
for i = 1:2
    for j = 1:2
        subplot(2,2,2*(i-1)+j); hold on
        plot([-1 1],[0 0],'k')
        plot([0 0],[-1 1],'k')
        % for i = 1:n
        for k = 1:n
            %         plot(real(B(i,j)),imag(B(i,j)),'.','MarkerSize',20,'Color',squeeze(map(i,j,:)))
            plot(real(B(i,j,k)),imag(B(i,j,k)),'.','MarkerSize',20,'Color',map(k,:))
        end
        % end
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
        % for i = 1:n
        for k = 1:n
            %         plot(real(B(i,j)),imag(B(i,j)),'.','MarkerSize',20,'Color',squeeze(map(i,j,:)))
            plot(real(F(i,j,k)),imag(F(i,j,k)),'.','MarkerSize',20,'Color',map(k,:))
        end
        % end
        if fixedHur
            title(['F (Hur = ' num2str(comp(idx)) ')'])
        else
            title(['F (Hud = ' num2str(comp(idx)) ')'])
        end
        axis equal
    end
end