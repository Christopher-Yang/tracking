% plot Bode plots for baseline and late learning data files

freqX = data{subj}{1}.sineParams.tX_freq{1};
freqY = data{subj}{1}.sineParams.tY_freq{1};

Nfreq = length(data{1}{1}.sineParams.tX_freq{1});
col = copper;
col = col(floor((size(col,1)/(Nfreq))*(1:Nfreq)),:);

names1 = {'xTarg_x','yTarg_x','xTarg_y','yTarg_y'};
names2 = {'xCurs_x','yCurs_x','xCurs_y','yCurs_y'};
clear Hur Hud Hur2 Hud2 B B2 F F2
for subj = 1:2
    for j = 1:4
        for i = 1:6
            a = data{subj}{i}.phasors;
            Hur(:,i,j,subj) = a.(names1{j}){1}.ratio;
            Hud(:,i,j,subj) = a.(names2{j}){2}.ratio;
            Hur2(:,i,j,subj) = reshape(permute([a.(names1{j}){3}.ratio a.(names1{j}){4}.ratio],[2 1]),[6 1]);
            Hud2(:,i,j,subj) = reshape(permute([a.(names2{j}){4}.ratio a.(names2{j}){3}.ratio],[2 1]),[6 1]);
        end
    end
end

Hur_gain = abs(Hur);
Hud_gain = abs(Hud);
Hur_phase = unwrap(angle(Hur));
Hud_phase = unwrap(angle(Hud));

for subj = 1:2
    count = 1;
    for k = 1:2
        for j = 1:2
            for i = 1:6
                a = data{subj}{i};
                B(:,i,count,subj) = a.B(j,k,:);
                B2(:,i,count,subj) = a.B2(j,k,:);
                F(:,i,count,subj) = a.F(j,k,:);
                F2(:,i,count,subj) = a.F2(j,k,:);
            end
            count = count + 1;
        end
    end
end

B_gain = abs(B);
F_gain = abs(F);
B_phase = unwrap(angle(B));
F_phase = unwrap(angle(F));

col = [1 0 0
       0 0 0];
gblocks = [1 3 2 4];
   
% figure(1); clf
% for j = 1:4
%     subplot(2,2,gblocks(j)); hold on
%     for i = 1:2
%         plot(freqX,Hur_gain(:,:,j,i),'Color',[col(i,:) 0.5])
%         plot(freqX,mean(Hur_gain(:,:,j,i),2),'Color',col(i,:),'LineWidth',2)
%     end
% end
% 
% figure(2); clf
% for j = 1:4
%     subplot(2,2,gblocks(j)); hold on
%     for i = 1:2
%         plot(freqX,B_gain(:,:,j,i),'Color',[col(i,:) 0.5])
%         plot(freqX,mean(B_gain(:,:,j,i),2),'Color',col(i,:),'LineWidth',2)
%     end
% end

figure(3); clf
for j = 1:4
    subplot(2,2,gblocks(j)); hold on
    for i = 1:2
        plot(freqX,Hur_phase(:,:,j,i),'Color',[col(i,:) 0.5])
        plot(freqX,mean(Hur_phase(:,:,j,i),2),'Color',col(i,:),'LineWidth',2)
    end
end

figure(4); clf
for j = 1:4
    subplot(2,2,gblocks(j)); hold on
    for i = 1:2
        plot(freqX,B_phase(:,:,j,i),'Color',[col(i,:) 0.5])
        plot(freqX,mean(B_phase(:,:,j,i),2),'Color',col(i,:),'LineWidth',2)
    end
end

figure(5); clf
for j = 1:4
    subplot(2,2,gblocks(j)); hold on
    for i = 1:2
        plot(freqX,F_phase(:,:,j,i),'Color',[col(i,:) 0.5])
        plot(freqX,mean(F_phase(:,:,j,i),2),'Color',col(i,:),'LineWidth',2)
    end
end