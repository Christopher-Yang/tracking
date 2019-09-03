function graph_coherence(data, gblocks, graph_name,output)
    
    rng(1);
    col = lines;
    col = col(1:7,:);
    
    a = {'x_x','y_y'};
    Nfreq = length(data{end}.ampX);
    Nblock = length(fieldnames(data{1}));
    SRboot = NaN(Nblock,Nfreq,1000);
    RRboot = NaN(Nblock,Nfreq,1000);
%     diffboot = NaN(Nblock,Nfreq,1000);
    for i = 1:length(a)
        dat = data{end}.(output).(a{i});
%         diff.(a{i}) = sqrt(dat.RRcohere) - dat.SRcohere;
        SRcohere = dat.SRcohere_full;
        RRcohere = dat.RRcohere_full;
%         diff_all = sqrt(dat.RRcohere_full) - dat.SRcohere_full;
        for j = 1:1000
            k = datasample(SRcohere,10,3);
            SRboot(:,:,j) = mean(k,3);
            k = datasample(sqrt(RRcohere),10,3);
            RRboot(:,:,j) = mean(k,3);
%             k = datasample(diff_all,10,3);
%             diffboot(:,:,j) = mean(k,3);
        end
        SRboot = sort(SRboot,3);
        RRboot = sort(RRboot,3);
%         diffboot = sort(diffboot,3);
        
        SRerr.(a{i}) = cat(3,SRboot(:,:,end-25),SRboot(:,:,26));
        SRerr.(a{i})(:,:,1) = SRerr.(a{i})(:,:,1) - dat.SRcohere;
        SRerr.(a{i})(:,:,2) = dat.SRcohere - SRerr.(a{i})(:,:,2);
        SRerr.(a{i}) = permute(SRerr.(a{i}), [3 2 1]);
        
        RRerr.(a{i}) = cat(3,RRboot(:,:,end-25),RRboot(:,:,26));
        RRerr.(a{i})(:,:,1) = RRerr.(a{i})(:,:,1) - sqrt(dat.RRcohere);
        RRerr.(a{i})(:,:,2) = sqrt(dat.RRcohere) - RRerr.(a{i})(:,:,2);
        RRerr.(a{i}) = permute(RRerr.(a{i}), [3 2 1]);
        
%         Derr.(a{i}) = cat(3,diffboot(:,:,end-25),diffboot(:,:,26));
%         Derr.(a{i})(:,:,1) = Derr.(a{i})(:,:,1) - diff.(a{i});
%         Derr.(a{i})(:,:,2) = diff.(a{i}) - Derr.(a{i})(:,:,2);
%         Derr.(a{i}) = permute(Derr.(a{i}), [3 2 1]);
    end
    
    dat = data{end}.(output);
    f_x = dat.x_x.freqs;
    f_y = dat.y_y.freqs;
    
    figure
    subplot(2,2,1)
    for j = 1:length(gblocks)
        s = shadedErrorBar(f_x,dat.x_x.SRcohere(gblocks(j),:),SRerr.x_x(:,:,gblocks(j)),'lineProps','-o');
        editErrorBar(s,col(j,:),1);
    end
    set(gca,'Xscale','log','box','off','TickDir','out')
    ylabel([output,' X movements'])
    axis([0.09 2.3 0 1])
    title('SR')
    yticks(0:0.2:1)
    
    subplot(2,2,3)
    for j = 1:length(gblocks)
        s = shadedErrorBar(f_y,dat.y_y.SRcohere(gblocks(j),:),SRerr.y_y(:,:,gblocks(j)),'lineProps','-o');
        editErrorBar(s,col(j,:),1);
    end
    set(gca,'Xscale','log','box','off','TickDir','out')
    xlabel('Frequency (Hz)')
    ylabel([output,' Y movements'])
    axis([0.09 2.3 0 1])
    yticks(0:0.2:1)
    
    subplot(2,2,2)
    for j = 1:length(gblocks)
        s = shadedErrorBar(f_x,sqrt(dat.x_x.RRcohere(gblocks(j),:)),RRerr.x_x(:,:,gblocks(j)),'lineProps','-o');
        editErrorBar(s,col(j,:),1);
    end
    set(gca,'Xscale','log','box','off','TickDir','out')
    axis([0.09 2.3 0 1])
    title('RR')
    yticks(0:0.2:1)
    legend(graph_name{gblocks},'Location','southwest')

    subplot(2,2,4)
    for j = 1:length(gblocks)
        s = shadedErrorBar(f_y,sqrt(dat.y_y.RRcohere(gblocks(j),:)),RRerr.y_y(:,:,gblocks(j)),'lineProps','-o');
        editErrorBar(s,col(j,:),1);
    end
    set(gca,'Xscale','log','box','off','TickDir','out')
    xlabel('Frequency (Hz)')
    axis([0.09 2.3 0 1])
    yticks(0:0.2:1)
    
%     subplot(2,3,3)
%     for j = 1:length(gblocks)
%         s = shadedErrorBar(f_x,diff.x_x(gblocks(j),:),Derr.x_x(:,:,gblocks(j)),'lineProps','-o');
%         editErrorBar(s,col(j,:),1);
%     end
%     set(gca,'Xscale','log','box','off','TickDir','out')
%     axis([0.09 2.3 0 1])
%     title('Unexplained')
%     
%     subplot(2,3,6)
%     for j = 1:length(gblocks)
%         s = shadedErrorBar(f_y,diff.y_y(gblocks(j),:),Derr.y_y(:,:,gblocks(j)),'lineProps','-o');
%         editErrorBar(s,col(j,:),1);
%     end
%     set(gca,'Xscale','log','box','off','TickDir','out')
%     axis([0.09 2.3 0 1])
%     legend(graph_name{gblocks})
end