function graph_traj(data)
% plots raw cursor and target trajectories in Figure 2B

% set variables for plotting
gblocks = [1:2 5:6];
block_name = {'baseline','pert1','pert2','pert3','pert4','post'}; 
graph_name = {'Baseline','Early','Late','Post'};
col = [166 124 82]./255;
Nblock = length(gblocks);

% plot rotation data
subj = 3;
trial = 3;
a = data.rot{subj};
figure(1); clf
for i = 1:Nblock % iterate over blocks
    subplot(2,4,i); hold on
    
    % plot target
    plot(a.(block_name{gblocks(i)}).target.x_pos(1:650,trial),a ...
        .(block_name{gblocks(i)}).target.y_pos(1:650,trial),'k');
    
    % plot cursor
    plot(a.(block_name{gblocks(i)}).cursor.x_pos(1:650,trial),a ...
        .(block_name{gblocks(i)}).cursor.y_pos(1:650,trial),'Color',col);
    
    pbaspect([1 1 1])
    axis([-.14 .14 -.14 .14])
    title(graph_name(i))
    if i == 4
        legend({'Target','Cursor'})
    elseif i == 1
        ylabel('Rotation: y position (m)')
    end
end

% plot mirror-reversal data
subj = 7;
trial = 5;
a = data.mir{subj};
for i = 1:Nblock % iterate over blocks
    subplot(2,4,i+4); hold on
    
    % plot 5 cm line as reference
    if i == 1
        plot([0 0.05], [-.1 -.1], 'k', 'LineWidth', 3)
    end
    
    % plot target
    plot(a.(block_name{gblocks(i)}).target.x_pos(1:650,trial),a ...
        .(block_name{gblocks(i)}).target.y_pos(1:650,trial),'k');

    % plot cursor
    plot(a.(block_name{gblocks(i)}).cursor.x_pos(1:650,trial),a ...
        .(block_name{gblocks(i)}).cursor.y_pos(1:650,trial),'Color',col);
    
    pbaspect([1 1 1])
    axis([-.14 .14 -.14 .14])
    if i == 1
        ylabel('Mirror Reversal: y position (m)')
    end
    xlabel('x position (m)')
end
end