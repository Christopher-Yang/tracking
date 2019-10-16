function graph_coherence(data, groups, block_name, gblocks, graph_name, output)
    
rng(1);
names1 = {'x_all','y_all'};
names2 = {'x','y'};
Nsubj = length(data.(groups{1}))-1;
Nfreq = length(data.(groups{1}){end}.ampX);
Nreps = length(data.(groups{1}){1}.(block_name{1}).MSE);
Nblock = length(block_name);
f_x = data.(groups{1}){end}.freqX;
f_y = data.(groups{1}){end}.freqY;
sorted_freqs = sort([f_x f_y]);

col = lines;
col = col(1:7,:);

for p = 1:length(groups)
    for i = 1:Nsubj
        for j = 1:Nblock
            targetX = mean(data.(groups{p}){i}.(block_name{j}).target.x_pos_all,2);
            targetY = mean(data.(groups{p}){i}.(block_name{j}).target.y_pos_all,2);
            outX = data.(groups{p}){i}.(block_name{j}).(output).x_pos_all;
            outY = data.(groups{p}){i}.(block_name{j}).(output).y_pos_all;
            
            N = size(targetX,1);
            for k = 1:Nreps
                coh = mscohere([outX(:,k) outY(:,k)],[targetX targetY],blackmanharris(round(N/5)),[],sorted_freqs,130.004,'mimo')';
                cohxx(k,:) = coh(1,:);
                cohyy(k,:) = coh(2,:);
            end
            SR.x_all(j,:,i,p) = mean(cohxx(:,1:2:end));
            SR.y_all(j,:,i,p) = mean(cohyy(:,2:2:end));
            
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
            RR.x_all(j,:,i,p) = mean(sqrt(cohx));
            RR.y_all(j,:,i,p) = mean(sqrt(cohy));
        end
    end
end

SR.x = squeeze(mean(SR.x_all,3));
SR.y = squeeze(mean(SR.y_all,3));
RR.x = squeeze(mean(RR.x_all,3));
RR.y = squeeze(mean(RR.y_all,3));

figure(1); clf
for j = 1:2
    subplot(2,2,j); hold on
    for i = 1:length(gblocks)
        plot(f_x,SR.x(gblocks(i),:,j),'-o','MarkerFaceColor',col(i,:),'MarkerEdgeColor','none')
    end
    if j == 1
        title('Rotation')
    else
        title('Mirror-Reversal')
    end
    ylabel([output,' SR_X coherence'])
    yticks(0:0.2:1)
    axis([0 2.3 0 1])
    
    subplot(2,2,j+2); hold on
    for i = 1:length(gblocks)
        plot(f_y,SR.y(gblocks(i),:,j),'-o','MarkerFaceColor',col(i,:),'MarkerEdgeColor','none')
    end
    xlabel('Frequency (Hz)')
    ylabel([output,' SR_Y coherence'])
    yticks(0:0.2:1)
    axis([0 2.3 0 1])
end
legend(graph_name(gblocks),'Location','southeast')
    
figure(2); clf
for j = 1:2
    subplot(2,2,j); hold on
    for i = 1:length(gblocks)
        plot(f_x,RR.x(gblocks(i),:,j),'-o','MarkerFaceColor',col(i,:),'MarkerEdgeColor','none')
    end
    if j == 1
        title('Rotation')
    else
        title('Mirror-Reversal')
    end
    ylabel([output,' RR_X coherence'])
    yticks(0:0.2:1)
    axis([0 2.3 0 1])
    
    subplot(2,2,j+2); hold on
    for i = 1:length(gblocks)
        plot(f_y,RR.y(gblocks(i),:,j),'-o','MarkerFaceColor',col(i,:),'MarkerEdgeColor','none')
    end
    ylabel([output,' RR_Y coherence'])
    xlabel('Frequency (Hz)')
    yticks(0:0.2:1)
    axis([0 2.3 0 1])
end
legend(graph_name(gblocks),'Location','southeast')
