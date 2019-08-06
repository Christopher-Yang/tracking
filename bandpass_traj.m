figure(1); clf
subj = 3;
block = 3;
time = data{subj}.(block_name{block}).time(1:end-13)/1000;
type = 'cursor';
lims = [0 0];

for i = 1:4
    if i == 1
        bp = [0.24 0.26];
        in = bandpass(data{subj}.(block_name{block}).target.x_pos,bp,130.004);
        out = bandpass(data{subj}.(block_name{block}).(type).x_pos,bp,130.004);
    elseif i == 2
        bp = [0.34 0.36];
        in = bandpass(data{subj}.(block_name{block}).target.y_pos,bp,130.004);
        out = bandpass(data{subj}.(block_name{block}).(type).x_pos,bp,130.004);
    elseif i == 3
        bp = [0.24 0.26];
        in = bandpass(data{subj}.(block_name{block}).target.x_pos,bp,130.004);
        out = bandpass(data{subj}.(block_name{block}).(type).y_pos,bp,130.004);
    else
        bp = [0.34 0.36];
        in = bandpass(data{subj}.(block_name{block}).target.y_pos,bp,130.004);
        out = bandpass(data{subj}.(block_name{block}).(type).y_pos,bp,130.004);
    end
    
    if lims(1) > min(in)+0.02*min(in) || lims(2) < max(in)+0.02*max(in)
        lims = [min(in)+0.02*min(in) max(in)+0.02*max(in)];
    end
    
    [pks,loc] = findpeaks(in);
    loc2 = [loc'; loc'];
    n = length(loc);
    x = [-ones(1,n); ones(1,n)];
    
    [pks,loc3] = findpeaks(out);
    loc4 = [loc3'; loc3'];
    n = length(loc3);
    x2 = [-ones(1,n); ones(1,n)];
    
    figure(1)
    subplot(2,2,i); hold on
    plot(time,in,'k')
    plot(time,out,'r')
    plot(time(loc2),x,'k')
    plot(time(loc4),x2,'r')
    ylim(lims)
    if i == 1
        ylabel('X_{Hand/Cursor}')
    elseif i == 3
        ylabel('Y_{Hand/Cursor}')
        xlabel('X_{Target}')
    elseif i == 4
        xlabel('Y_{Target}')
    end
end

for i = 1:4
    subplot(2,2,i)
    ylim(lims)
end