function graph_coherence(data, gblocks, graph_name)
    
    rng(1);
    col = lines;
    col = col(1:7,:);
    a = {'x_x','y_y'};
    Nfreq = length(data{end}.x_x.freqs);
    Nblock = size(data{end}.x_x.amplitude,1);
%     SRboot = NaN(Nblock,Nfreq*2,1000);
    SRboot = NaN(Nblock,Nfreq,1000);
    RRboot = NaN(Nblock,Nfreq,1000);
    diffboot = NaN(Nblock,Nfreq,1000);
    for i = 1:length(a)
        diff.(a{i}) = sqrt(data{end}.(a{i}).RRcohere) - data{end}.(a{i}).SRcohere;
        SRcohere = data{end}.(a{i}).SRcohere_full;
        RRcohere = data{end}.(a{i}).RRcohere_full;
        diff_all = sqrt(data{end}.(a{i}).RRcohere_full) - data{end}.(a{i}).SRcohere_full;
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
        
        SRerr.(a{i}) = cat(3,SRboot(:,:,end-25),SRboot(:,:,26));
        SRerr.(a{i})(:,:,1) = SRerr.(a{i})(:,:,1) - data{end}.(a{i}).SRcohere;
        SRerr.(a{i})(:,:,2) = data{end}.(a{i}).SRcohere - SRerr.(a{i})(:,:,2);
        SRerr.(a{i}) = permute(SRerr.(a{i}), [3 2 1]);
        
        RRerr.(a{i}) = cat(3,RRboot(:,:,end-25),RRboot(:,:,26));
        RRerr.(a{i})(:,:,1) = RRerr.(a{i})(:,:,1) - sqrt(data{end}.(a{i}).RRcohere);
        RRerr.(a{i})(:,:,2) = sqrt(data{end}.(a{i}).RRcohere) - RRerr.(a{i})(:,:,2);
        RRerr.(a{i}) = permute(RRerr.(a{i}), [3 2 1]);
        
        Derr.(a{i}) = cat(3,diffboot(:,:,end-25),diffboot(:,:,26));
        Derr.(a{i})(:,:,1) = Derr.(a{i})(:,:,1) - diff.(a{i});
        Derr.(a{i})(:,:,2) = diff.(a{i}) - Derr.(a{i})(:,:,2);
        Derr.(a{i}) = permute(Derr.(a{i}), [3 2 1]);
    end
    
    f_x = data{end}.x_x.freqs;
    f_y = data{end}.y_y.freqs;
    freqs_all = sort([f_x f_y]);
    names = {'Rotation','Mirror Reversal'};
    
    figure('Name','Linearity','NumberTitle','off');
    for i = 1:length(groups)
        subplot(2,2,i)
        for j = 1:length(gblocks)
            s = shadedErrorBar(f_x,data{end}.x_x.SRcohere(gblocks(j),:),SRerr{i}.x_x(:,:,gblocks(j)),'lineProps','-o');
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
            s = shadedErrorBar(f_y,data{end}.y_y.SRcohere(gblocks(j),:),SRerr{i}.x_x(:,:,gblocks(j)),'lineProps','-o');
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
            s = shadedErrorBar(f_x,sqrt(data{end}.x_x.RRcohere(gblocks(j),:)),RRerr{i}.x_x(:,:,gblocks(j)),'lineProps','-o');
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
            s = shadedErrorBar(f_y,sqrt(data{end}.y_y.RRcohere(gblocks(j),:)),RRerr{i}.x_x(:,:,gblocks(j)),'lineProps','-o');
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