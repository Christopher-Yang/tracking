col = lines;
col = col(1:7,:);amp = NaN(2,2,420);
amp2 = NaN(2,2,420);
z = NaN(420,1);
y = NaN(420,1);
z_boot = NaN(6,7,1000);
y_boot = NaN(6,7,1000);
z_err = NaN(6,7,2);
y_err = NaN(6,7,2);

amp(1,1,:) = reshape(data.rot.avg.x_x.all_amp,[1 1 420]);
amp(2,2,:) = reshape(data.rot.avg.y_y.all_amp,[1 1 420]);
amp(1,2,:) = reshape(data.rot.avg.x_y.all_amp,[1 1 420]);
amp(2,1,:) = reshape(data.rot.avg.y_x.all_amp,[1 1 420]);

amp2(1,1,:) = reshape(data.rot_i.avg.x_x.all_amp,[1 1 420]);
amp2(2,2,:) = reshape(data.rot_i.avg.y_y.all_amp,[1 1 420]);
amp2(1,2,:) = reshape(data.rot_i.avg.x_y.all_amp,[1 1 420]);
amp2(2,1,:) = reshape(data.rot_i.avg.y_x.all_amp,[1 1 420]);

for i = 1:420
    z(i,1) = det(amp(:,:,i));
    y(i,1) = det(amp2(:,:,i));
end

z2 = reshape(z,[6,7,10]);
y2 = reshape(y,[6 7,10]);
z_avg = mean(z2,3);
y_avg = mean(y2,3);

for i = 1:1000
    p = datasample(z2,10,3);
    s = datasample(y2,10,3);
    z_boot(:,:,i) = mean(p,3);
    y_boot(:,:,i) = mean(s,3);
end

z_boot = sort(z_boot,3);
y_boot = sort(y_boot,3);

z_err(:,:,2) = z_avg - z_boot(:,:,26);
z_err(:,:,1) = z_boot(:,:,end-25) - z_avg;
y_err(:,:,2) = y_avg - y_boot(:,:,26);
y_err(:,:,1) = y_boot(:,:,end-25) - y_avg;

z_err = permute(z_err,[3 2 1]);
y_err = permute(y_err,[3 2 1]);

% f_x = data.rot.avg.x_x.freqs;
% f_y = data.rot.avg.y_y.freqs;

z_avg = [z_avg(1,:); z_avg(2:5,:); z_avg(6,:)];
y_avg = [y_avg(1,:); y_avg(2:5,:); y_avg(6,:)];

gblocks = [1:2 5:6];
figure;
subplot(1,2,1);
% semilogx(f_x,z2(:,[1 2 5 6]),'LineWidth',3);
for i = 1:length(gblocks)
    zz = shadedErrorBar(1:7,z_avg(gblocks(i),:),z_err(:,:,gblocks(i)),'lineProps','-ro');
    editErrorBar(zz,col(i,:),3);
end
title('Rotation');
ylabel('Determinant');
xticks([1 7]);
xticklabels({'Lowest Frequency','Highest Frequency'});

subplot(1,2,2);
for i = 1:length(gblocks)
    zz = shadedErrorBar(1:7,y_avg(gblocks(i),:),y_err(:,:,gblocks(i)),'lineProps','-ro');
    editErrorBar(zz,col(i,:),3);
end
title('Mirror Reversal');
ylabel('Determinant');
xticks([1 7]);
xticklabels({'Lowest Frequency','Highest Frequency'});
legend({'Baseline','Naive','Max Training','Aftereffects'});