interval = pi/16;
angle = pi:interval:2*pi-interval;
angle = fliplr(angle);
realPart = cos(angle);
imaginaryPart = sin(angle);
comp = permute(realPart + imaginaryPart*1j,[2 1]);
template = [real(comp)'; imag(comp)'];
amp = linspace(0.8,0.1,length(comp))';
comp = 0.7.*comp;

clear Hud Hur B F
Hud(1,1,:) = comp(end-15:end);
Hud(2,2,:) = comp(end-15:end);
Hur(1,1,:) = comp(1:16);
Hur(2,2,:) = comp(1:16);
M = eye(2);

n = size(Hur,3);

% disturbance after M
for i = 1:size(Hud,3)
    B(:,:,i) = -Hud(:,:,i)*inv(M*Hud(:,:,i) + eye(2));
    F(:,:,i) = Hur(:,:,i) + B(:,:,i)*(M*Hur(:,:,i) - eye(2));
end

% set the visuomotor perturbation
M2 = rotz(-90);
M2 = M2(1:2,1:2);

% set the degree to which the controllers rotate
R_B = rotz(40);
R_B = R_B(1:2,1:2);

R_F = rotz(40);
R_F = R_F(1:2,1:2);

clear Hur_rot Hud_rot
for i = 1:n
    B_rot = R_B*B(:,:,i);
    F_rot = R_F*F(:,:,i);
    
    Hur_rot(:,:,i) = inv(eye(2) + B_rot*M2)*(F_rot + B_rot);
    Hud_rot(:,:,i) = -inv(eye(2) + B_rot*M2)*B_rot;
end

clear Hur_gain Hud_gain
for i = 1:2
    for j = 1:2
        % compute gain matrix for Hur
        phasors = squeeze(Hur_rot(i,j,:));
        phasors = [real(phasors) imag(phasors)];
        
        clear gain
        for k = 1:n
            gain(k) = phasors(k,:)*template(:,k);
        end
        Hur_gain(i,j,:) = gain;
        
        % compute gain matrix for Hud
        phasors = squeeze(Hud_rot(i,j,:));
        phasors = [real(phasors) imag(phasors)];
        
        clear gain
        for k = 1:n
            gain(k) = phasors(k,:)*template(:,k);
        end
        Hud_gain(i,j,:) = gain;
    end
end

% for plotting heatmaps
col1 = [1 0 0];
col2 = [1 1 1];
Nstep = 100;
map1 = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2),Nstep)', linspace(col1(3),col2(3),Nstep)'];

col1 = [1 1 1];
col2 = [0 0 1];
map2 = [linspace(col1(1),col2(1),Nstep)', linspace(col1(2),col2(2),Nstep)', linspace(col1(3),col2(3),Nstep)'];

map = [map1; map2];
clims = [-1 1];

figure(1); clf
for i = 1:n
    subplot(2,n,i)
    imagesc(Hur_gain(:,:,i),clims)
    colormap(map)
    pbaspect([1 1 1])
    set(gca,'Xtick',[],'Ytick',[])
    if i == 1
        ylabel('Hur')
    end
    
    subplot(2,n,i+n)
    imagesc(Hud_gain(:,:,i),clims)
    colormap(map)
    pbaspect([1 1 1])
    set(gca,'Xtick',[],'Ytick',[])
    if i == 1
        ylabel('Hud')
    end
end

nReal = numel(realPart);
nImag = numel(imaginaryPart);
col1 = [1 0 0];
col2 = [0 1 0];
col3 = [0 0 1];
col4 = [1 1 0];

map2 = [linspace(col1(1),col4(1),n)' linspace(col1(2),col4(2),n)' linspace(col1(3),col4(3),n)'];

% figure(2); clf
% for i = 1:2
%     for k = 1:2
%         subplot(2,2,2*(i-1)+k); hold on
%         plot([-1 1],[0 0],'k')
%         plot([0 0],[-1 1],'k')
%             for j = 1:n
%                 plot(real(Hud_rot(i,k,j)),imag(Hud_rot(i,k,j)),'.','MarkerSize',20,'Color',map2(j,:))
%             end
%         title('H_{ud}')
%         axis equal
%     end
% end
% 
% figure(3); clf
% for i = 1:2
%     for k = 1:2
%         subplot(2,2,2*(i-1)+k); hold on
%         plot([-1 1],[0 0],'k')
%         plot([0 0],[-1 1],'k')
%             for j = 1:n
%                 plot(real(Hur_rot(i,k,j)),imag(Hur_rot(i,k,j)),'.','MarkerSize',20,'Color',map2(j,:))
%             end
%         title('H_{ur}')
%         axis equal
%     end
% end

%%
% figure(8); clf
% for i = 1:2
%     for j = 1:2
%         subplot(2,2,2*(i-1)+j); hold on
%         plot([-1 1],[0 0],'k')
%         plot([0 0],[-1 1],'k')
%         % for i = 1:n
%         for k = 1:n
%             %         plot(real(B(i,j)),imag(B(i,j)),'.','MarkerSize',20,'Color',squeeze(map(i,j,:)))
%             plot(real(B_rot(i,j,k)),imag(B_rot(i,j,k)),'.','MarkerSize',20,'Color',map(k,:))
%         end
%         % end
%         title('B')
%         axis equal
%     end
% end
% 
% figure(9); clf
% for i = 1:2
%     for j = 1:2
%         subplot(2,2,2*(i-1)+j); hold on
%         plot([-1 1],[0 0],'k')
%         plot([0 0],[-1 1],'k')
%         % for i = 1:n
%         for k = 1:n
%             %         plot(real(B(i,j)),imag(B(i,j)),'.','MarkerSize',20,'Color',squeeze(map(i,j,:)))
%             plot(real(F_rot(i,j,k)),imag(F_rot(i,j,k)),'.','MarkerSize',20,'Color',map(k,:))
%         end
%         % end
%         if fixedHur
%             title(['F (Hur = ' num2str(comp(idx)) ')'])
%         else
%             title(['F (Hud = ' num2str(comp(idx)) ')'])
%         end
%         axis equal
%     end
% end

%%

nReal = numel(realPart);
nImag = numel(imaginaryPart);
col1 = [1 0 0];
col2 = [0 1 0];
col3 = [0 0 1];
col4 = [1 1 0];

map = [linspace(col1(1),col4(1),n)' linspace(col1(2),col4(2),n)' linspace(col1(3),col4(3),n)'];

ax1 = 1;
ax2 = 1;

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
        title('F')
        axis equal
    end
end