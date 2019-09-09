function graph_coherence(data, block_name, graph_name, output, subj)
    
rng(1);
names1 = {'x_all','y_all'};
names2 = {'x','y'};
Nsubj = length(subj);
Nfreq = length(data{end}.ampX);
Nreps = length(data{1}.(block_name{1}).MSE);
Nblock = length(block_name);
f_x = data{end}.cursor.x_x.freqs;
f_y = data{end}.cursor.y_y.freqs;
sorted_freqs = sort([f_x f_y]);

col = copper;
n = round(64.*((1:Nfreq-1)/(Nfreq-1)));
col = col([1 n],:);

for i = 1:Nsubj
    for j = 1:Nblock
        targetX = mean(data{subj(i)}.(block_name{j}).target.x_pos_all,2);
        targetY = mean(data{subj(i)}.(block_name{j}).target.y_pos_all,2);
        outX = data{subj(i)}.(block_name{j}).(output).x_pos_all;
        outY = data{subj(i)}.(block_name{j}).(output).y_pos_all;
        
        N = size(targetX,1);
        for k = 1:Nreps
            coh = mscohere([outX(:,k) outY(:,k)],[targetX targetY],blackmanharris(round(N/5)),[],sorted_freqs,130.004,'mimo')';
            cohxx(k,:) = coh(1,:);
            cohyy(k,:) = coh(2,:);
        end
        SR.x_all(j,:,i) = mean(cohxx(:,1:2:end));
        SR.y_all(j,:,i) = mean(cohyy(:,2:2:end));
        
        FLAG = 1;
        k1 = 1;
        k2 = Nreps;
        k3 = 1;
        while FLAG
            cohx(k3,:) = mscohere(outX(:,k1),outX(:,k2),[],[],f_x,130.004);
            cohy(k3,:) = mscohere(outY(:,k1),outY(:,k2),[],[],f_y,130.004);
            k1 = k1 + 1;
            if k1 == k2
                k1 = 1;
                k2 = k2 - 1;
            end
            if k2 == 1
                FLAG = 0;
            end
            k3 = k3 + 1;
        end
        RR.x_all(j,:,i) = mean(cohx);
        RR.y_all(j,:,i) = mean(cohy);
    end
end

SR.x = mean(SR.x_all,3);
SR.y = mean(SR.y_all,3);
RR.x = mean(RR.x_all,3);
RR.y = mean(RR.y_all,3);

% SRboot = NaN(Nblock,Nfreq,1000);
% RRboot = NaN(Nblock,Nfreq,1000);
% for i = 1:length(names1)
%     SRcohere = SR.(names1{i});
%     RRcohere = RR.(names1{i});
%     for j = 1:1000
%         k = datasample(SRcohere,10,3);
%         SRboot(:,:,j) = mean(k,3);
%         k = datasample(sqrt(RRcohere),10,3);
%         RRboot(:,:,j) = mean(k,3);
%     end
%     SRboot = sort(SRboot,3);
%     RRboot = sort(RRboot,3);
%     
%     SRerr.(names2{i}) = cat(3,SRboot(:,:,end-25),SRboot(:,:,26));
%     SRerr.(names2{i})(:,:,1) = SRerr.(names2{i})(:,:,1) - SR.(names2{i});
%     SRerr.(names2{i})(:,:,2) = SR.(names2{i}) - SRerr.(names2{i})(:,:,2);
%     SRerr.(names2{i}) = permute(SRerr.(names2{i}), [3 2 1]);
%     
%     RRerr.(names2{i}) = cat(3,RRboot(:,:,end-25),RRboot(:,:,26));
%     RRerr.(names2{i})(:,:,1) = RRerr.(names2{i})(:,:,1) - sqrt(RR.(names2{i}));
%     RRerr.(names2{i})(:,:,2) = sqrt(RR.(names2{i})) - RRerr.(names2{i})(:,:,2);
%     RRerr.(names2{i}) = permute(RRerr.(names2{i}), [3 2 1]);
% end

% figure
% subplot(2,2,1); hold on
% % for j = 1:length(gblocks)
% %     s = shadedErrorBar(f_x,SR.x(j,:),SRerr.x(:,:,j),'lineProps','-o');
% %     editErrorBar(s,col(j,:),1);
% % end
% for i = 1:Nblock
%     plot(f_x, SR.x(i,:),'-o','Color',col(i,:),'LineWidth',1,'MarkerFaceColor',col(i,:),'MarkerEdgeColor','none')
%     plot(f_x,squeeze(SR.x_all(i,:,:)),'Color',[col(i,:) 0.4],'HandleVisibility','off')
% end
% set(gca,'Xscale','log','box','off','TickDir','out')
% ylabel([output,' X movements'])
% axis([0.09 2.3 0 1])
% title('SR')
% yticks(0:0.2:1)
% legend(graph_name{gblocks},'Location','southwest')
% 
% subplot(2,2,3); hold on
% % for j = 1:length(gblocks)
% %     s = shadedErrorBar(f_y,SR.y(j,:),SRerr.y(:,:,j),'lineProps','-o');
% %     editErrorBar(s,col(j,:),1);
% % end
% for i = 1:Nblock
%     plot(f_y, SR.y(i,:),'-o','Color',col(i,:),'LineWidth',1,'MarkerFaceColor',col(i,:),'MarkerEdgeColor','none')
%     plot(f_y,squeeze(SR.y_all(i,:,:)),'Color',[col(i,:) 0.4])
% end
% set(gca,'Xscale','log','box','off','TickDir','out')
% xlabel('Frequency (Hz)')
% ylabel([output,' Y movements'])
% axis([0.09 2.3 0 1])
% yticks(0:0.2:1)
% 
% subplot(2,2,2); hold on
% % for j = 1:length(gblocks)
% %     s = shadedErrorBar(f_x,sqrt(RR.x(j,:)),RRerr.x(:,:,j),'lineProps','-o');
% %     editErrorBar(s,col(j,:),1);
% % end
% for i = 1:Nblock
%     plot(f_x, RR.x(i,:),'-o','Color',col(i,:),'LineWidth',1,'MarkerFaceColor',col(i,:),'MarkerEdgeColor','none')
%     plot(f_x,squeeze(RR.x_all(i,:,:)),'Color',[col(i,:) 0.4])
% end
% set(gca,'Xscale','log','box','off','TickDir','out')
% axis([0.09 2.3 0 1])
% title('RR')
% yticks(0:0.2:1)
% 
% subplot(2,2,4); hold on
% % for j = 1:length(gblocks)
% %     s = shadedErrorBar(f_y,sqrt(RR.y(j,:)),RRerr.y(:,:,j),'lineProps','-o');
% %     editErrorBar(s,col(j,:),1);
% % end
% for i = 1:Nblock
%     plot(f_y, RR.x(i,:),'-o','Color',col(i,:),'LineWidth',1,'MarkerFaceColor',col(i,:),'MarkerEdgeColor','none')
%     plot(f_y,squeeze(RR.y_all(i,:,:)),'Color',[col(i,:) 0.4])
% end
% set(gca,'Xscale','log','box','off','TickDir','out')
% xlabel('Frequency (Hz)')
% axis([0.09 2.3 0 1])
% yticks(0:0.2:1)

for i = 1:Nfreq
    leg{i} = ['Freq ',num2str(i)];
end
tick = find(contains(graph_name,'('));

figure
subplot(2,2,1); hold on
for i = 1:Nfreq
    plot(1:Nblock,SR.x(:,i),'-o','Color',col(i,:),'LineWidth',1,'MarkerFaceColor',col(i,:),'MarkerEdgeColor','none')
end
title('SR')
ylabel([output,' X coherence'])
xticks(tick)
xticklabels(graph_name(tick))
yticks(0:0.2:1)
ylim([0 1])
legend(leg,'Location','southeast')

subplot(2,2,3); hold on
for i = 1:Nfreq
    plot(1:Nblock,SR.y(:,i),'-o','Color',col(i,:),'LineWidth',1,'MarkerFaceColor',col(i,:),'MarkerEdgeColor','none')
end
xlabel('Block')
ylabel([output,' Y coherence'])
xticks(tick)
xticklabels(graph_name(tick))
yticks(0:0.2:1)
ylim([0 1])

subplot(2,2,2); hold on
for i = 1:Nfreq
    plot(1:Nblock,RR.x(:,i),'-o','Color',col(i,:),'LineWidth',1,'MarkerFaceColor',col(i,:),'MarkerEdgeColor','none')
end
title('RR')
xticks(tick)
xticklabels(graph_name(tick))
yticks(0:0.2:1)
ylim([0 1])

subplot(2,2,4); hold on
for i = 1:Nfreq
    plot(1:Nblock,RR.y(:,i),'-o','Color',col(i,:),'LineWidth',1,'MarkerFaceColor',col(i,:),'MarkerEdgeColor','none')
end
xlabel('Block')
xticks(tick)
xticklabels(graph_name(tick))
yticks(0:0.2:1)
ylim([0 1])
end