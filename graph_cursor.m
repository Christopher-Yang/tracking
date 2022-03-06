% graph_cursor plots example cursor and target trajectories from different
% blocks of the experiment

function graph_cursor(data)

time = 1:650; % 5 seconds of data to plot
names = {'Baseline','Early','Late (2-day)','Late (5-day)','Late (10-day)'};

figure(1); clf
for i = 1:5
    
    % select subject and block to plot
    switch i
        case 1
            a = data.day2{1}.B1_baseline;
        case 2
            a = data.day2{1}.B3;
        case 3
            a = data.day2{1}.B9;
        case 4
            a = data.day5{2}.B24;
        case 5
            a = data.day10{4}.B49;
    end
    
    subplot(1,5,i); hold on
    plot(a.target.x_pos(time,1), a.target.y_pos(time,1),'r');
    plot(a.cursor.x_pos(time,1), a.cursor.y_pos(time,1),'k');
    if i == 1
        plot([0 0.1],[-0.15 -0.15],'k','LineWidth',2)
    end
    axis([-0.15 0.15 -0.2 0.11])
    axis square
    title(names{i})
    xticks([])
    yticks([])
end

% save figure for Illustrator
% print('C:/Users/Chris/Documents/Papers/habit/figure_drafts/cursor','-dpdf','-painters')

end
