function graph_coherence(data, groups, gblocks, graph_name)
    
    rng(1);
    col = lines;
    col = col(1:7,:);
    a = {'x_x','y_y'};
    Nfreq = length(data.(groups{1}).avg.x_x.freqs);
    Nblock = size(data.(groups{1}).avg.x_x.amplitude,1);
%     SRboot = NaN(Nblock,Nfreq*2,1000);
    SRboot = NaN(Nblock,Nfreq,1000);
    RRboot = NaN(Nblock,Nfreq,1000);
    diffboot = NaN(Nblock,Nfreq,1000);
    for z = 1:length(groups)
        for i = 1:length(a)
            diff{z}.(a{i}) = sqrt(data.(groups{z}).avg.(a{i}).RRcohere) - data.(groups{z}).avg.(a{i}).SRcohere;
            SRcohere = data.(groups{z}).avg.(a{i}).SRcohere_full;
            RRcohere = data.(groups{z}).avg.(a{i}).RRcohere_full;
            diff_all = sqrt(data.(groups{z}).avg.(a{i}).RRcohere_full) - data.(groups{z}).avg.(a{i}).SRcohere_full;
            for j = 1:1000
                k = datasample(SRcohere,10,3);
                SRboot(:,:,j) = mean(k,3);
                k = datasample(sqrt(RRcohere),10,3);
                RRboot(:,:,j) = mean(k,3);
                k = datasample(diff_all,10,3);
                diffboot(:,:,j) = mean(k,3);
            end
            SRboot = sort(SRboot,3);
            RRboot = sort(RRboot,3);
            diffboot = sort(diffboot,3);
            
            SRerr{z}.(a{i}) = cat(3,SRboot(:,:,end-25),SRboot(:,:,26));
            SRerr{z}.(a{i})(:,:,1) = SRerr{z}.(a{i})(:,:,1) - data.(groups{z}).avg.(a{i}).SRcohere;
            SRerr{z}.(a{i})(:,:,2) = data.(groups{z}).avg.(a{i}).SRcohere - SRerr{z}.(a{i})(:,:,2);
            SRerr{z}.(a{i}) = permute(SRerr{z}.(a{i}), [3 2 1]);
            
            RRerr{z}.(a{i}) = cat(3,RRboot(:,:,end-25),RRboot(:,:,26));
            RRerr{z}.(a{i})(:,:,1) = RRerr{z}.(a{i})(:,:,1) - sqrt(data.(groups{z}).avg.(a{i}).RRcohere);
            RRerr{z}.(a{i})(:,:,2) = sqrt(data.(groups{z}).avg.(a{i}).RRcohere) - RRerr{z}.(a{i})(:,:,2);
            RRerr{z}.(a{i}) = permute(RRerr{z}.(a{i}), [3 2 1]);
            
            Derr{z}.(a{i}) = cat(3,diffboot(:,:,end-25),diffboot(:,:,26));
            Derr{z}.(a{i})(:,:,1) = Derr{z}.(a{i})(:,:,1) - diff{z}.(a{i});
            Derr{z}.(a{i})(:,:,2) = diff{z}.(a{i}) - Derr{z}.(a{i})(:,:,2);
            Derr{z}.(a{i}) = permute(Derr{z}.(a{i}), [3 2 1]);
        end
    end
    
    f_x = data.(groups{1}).avg.x_x.freqs;
    f_y = data.(groups{1}).avg.y_y.freqs;
    freqs_all = sort([f_x f_y]);
    names = {'Rotation','Mirror Reversal'};
    
    figure('Name','Linearity','NumberTitle','off');
    for i = 1:length(groups)
        subplot(2,2,i)
        for j = 1:length(gblocks)
            s = shadedErrorBar(f_x,data.(groups{i}).avg.x_x.SRcohere(gblocks(j),:),SRerr{i}.x_x(:,:,gblocks(j)),'lineProps','-o');
            editErrorBar(s,col(j,:),1);
        end
        set(gca,'Xscale','log','box','off','TickDir','out')
        if i == 1
            ylabel('SR_X Coherence')
        end
        axis([0.09 2.3 0.25 1])
        title(names{i})
        yticks(0:0.25:1)
        
        subplot(2,2,i+2)
        for j = 1:length(gblocks)
            s = shadedErrorBar(f_y,data.(groups{i}).avg.y_y.SRcohere(gblocks(j),:),SRerr{i}.x_x(:,:,gblocks(j)),'lineProps','-o');
            editErrorBar(s,col(j,:),1);
        end
        set(gca,'Xscale','log','box','off','TickDir','out')
        xlabel('Frequency (Hz)')
        if i == 1
            ylabel('SR_Y Coherence')
        end
        axis([0.09 2.3 0.25 1])
        yticks(0:0.25:1)
    end
    legend(graph_name{gblocks},'Location','southeast')
    
    figure('Name','Determinism','NumberTitle','off');
    for i = 1:length(groups)
%         RRerr{i}.x_x = RRerr{i}.x_x([2 1],:,:); % use for 1 - sqrt(RR)
        subplot(2,2,i)
        for j = 1:length(gblocks)
            s = shadedErrorBar(f_x,sqrt(data.(groups{i}).avg.x_x.RRcohere(gblocks(j),:)),RRerr{i}.x_x(:,:,gblocks(j)),'lineProps','-o');
            editErrorBar(s,col(j,:),1);
        end
        set(gca,'Xscale','log','box','off','TickDir','out','Ytick',0.25:0.25:1)
        if i == 1
            ylabel('$\sqrt{RR_X\: Coherence}$','Interpreter','latex')
            legend(graph_name{gblocks},'Location','southwest')
        end
        axis([0.09 2.3 0.25 1])
        title(names{i})
        pbaspect([1 1 1])
        
        RRerr{i}.y_y = RRerr{i}.y_y(:,:,[2 1]);
        subplot(2,2,i+2)
        for j = 1:length(gblocks)
            s = shadedErrorBar(f_y,sqrt(data.(groups{i}).avg.y_y.RRcohere(gblocks(j),:)),RRerr{i}.x_x(:,:,gblocks(j)),'lineProps','-o');
            editErrorBar(s,col(j,:),1);
        end
        set(gca,'Xscale','log','box','off','TickDir','out','Ytick',0.25:0.25:1)
        xlabel('Frequency (Hz)')
        if i == 1
            ylabel('$\sqrt{RR_Y\: Coherence}$','Interpreter','latex')
        end
        axis([0.09 2.3 0.25 1])
        pbaspect([1 1 1])
    end
    
    figure('Name','Nonlinearity','NumberTitle','off')
    for i = 1:length(groups)
        subplot(2,2,i)
        for j = 1:length(gblocks)
            s = shadedErrorBar(f_x,diff{i}.x_x(gblocks(j),:),Derr{i}.x_x(:,:,gblocks(j)),'lineProps','-o');
            editErrorBar(s,col(j,:),1);
        end
        set(gca,'Xscale','log','box','off','TickDir','out')
        if i == 1
            ylabel('\surd{RR_X} - SR_X Coherence')
        end
        axis([0.09 2.3 0 1])
        title(names{i})
        
        subplot(2,2,i+2)
        for j = 1:length(gblocks)
            s = shadedErrorBar(f_x,diff{i}.y_y(gblocks(j),:),Derr{i}.y_y(:,:,gblocks(j)),'lineProps','-o');
            editErrorBar(s,col(j,:),1);
        end
        set(gca,'Xscale','log','box','off','TickDir','out')
        if i == 1
            ylabel(['\surd{RR_Y} ' char(8211) ' SR_Y Coherence'])
        end
        xlabel('Frequency (Hz)')
        axis([0.09 2.3 0 1])
    end
    legend(graph_name{gblocks})
end