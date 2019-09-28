rng(1);
a = {'x_x','y_y'};
Nfreq = length(data.(groups{1}){end}.freqX);
Nblock = size(data.(groups{1}){end}.cursor.x_x.gain,1);
SRboot = NaN(Nblock,Nfreq,1000);
for z = 1:length(groups)
    for i = 1:length(a)
        SRcohere = data.(groups{z}){end}.(a{i}).SRcohere_full;
        for j = 1:1000
            k = datasample(SRcohere,10,3);
            SRboot(:,:,j) = mean(k,3);
        end
        SRboot = sort(SRboot,3);
        
        SRerr{z}.(a{i}) = cat(3,SRboot(:,:,end-25),SRboot(:,:,26));
        SRerr{z}.(a{i})(:,:,1) = SRerr{z}.(a{i})(:,:,1) - data.(groups{z}).avg.(a{i}).SRcohere;
        SRerr{z}.(a{i})(:,:,2) = data.(groups{z}).avg.(a{i}).SRcohere - SRerr{z}.(a{i})(:,:,2);
        SRerr{z}.(a{i}) = permute(SRerr{z}.(a{i}), [3 2 1]);
    end
end

f_x = data.(groups{1}).avg.x_x.freqs;
f_y = data.(groups{1}).avg.y_y.freqs;
freqs_all = sort([f_x f_y]);
names = {'Rotation','Mirror Reversal'};

col = [255 69 0
       65 105 225]/255;

figure; hold on
s = shadedErrorBar(f_x,data.rot.avg.x_x.SRcohere(1,:),SRerr{1}.x_x(:,:,1),'lineProps','-o');
editErrorBar(s,col(1,:),1);
s = shadedErrorBar(f_x,data.rot_i.avg.x_x.SRcohere(1,:),SRerr{2}.x_x(:,:,1),'lineProps','-o');
editErrorBar(s,col(2,:),1);
s = shadedErrorBar(f_y,data.rot.avg.y_y.SRcohere(1,:),SRerr{1}.x_x(:,:,1),'lineProps','--o');
editErrorBar(s,col(1,:),1);
s = shadedErrorBar(f_y,data.rot_i.avg.y_y.SRcohere(1,:),SRerr{2}.x_x(:,:,1),'lineProps','--o');
editErrorBar(s,col(2,:),1);
set(gca,'box','off','TickDir','out')
axis square
ylabel('Coherence')
axis([0 2.2 0.7 1])
yticks(0.7:0.1:1)
legend({'VMR_X','MR_X','VMR_Y','MR_Y'},'Location','southwest')