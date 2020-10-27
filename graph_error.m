freqX = data{1}{1}.sineParams.tX_freq{1};
freqY = data{1}{1}.sineParams.tY_freq{1};
Nblock = 6;
Nfreq = 6;

Ncols = 6;
col = copper;
col = col(floor((size(col,1)/(Ncols))*(1:Ncols)),:);

clear B B2 F F2
for i = 1:Nblock
    a = data{1}{i};
    count = 1;
    for j = 1:2
        for k = 1:2
            B(:,i,count) = a.B(j,k,:);
            B2(:,i,count) = a.B2(j,k,:);
            F(:,i,count) = a.F(j,k,:);
            F2(:,i,count) = a.F2(j,k,:);
            count = count + 1;
        end
    end
end

clear Berror Ferror
for i = 1:4
    Berror(:,:,i) = abs(B(:,1,i) - B(:,2:end,i));
    Ferror(:,:,i) = abs(F(:,1,i) - B(:,2:end,i));
end

figure(1); clf
for j = 1:4
    subplot(2,2,j); hold on
    for i = [1:3 6]
        plot(permute(Berror(i,:,j),[2 1]),'Color',col(i,:),'LineWidth',2)
    end
    xticks(1:5)
    xlabel('Block')
    ylabel('Error')
end

figure(2); clf
for j = 1:4
    subplot(2,2,j); hold on
    for i = [1:3 6]
        plot(permute(Ferror(i,:,j),[2 1]),'Color',col(i,:),'LineWidth',2)
    end
    xticks(1:5)
    xlabel('Block')
    ylabel('Error')
end
