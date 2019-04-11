subj = 2;
trial = 4;
a = d.(subj_rot_i{subj});
q = [1:2 5:6];

figure
for i = 1:4
    subplot(2,2,i)
    plot(a.(block_name{q(i)}).traj(:,1,trial),a.(block_name{q(i)}).traj(:,2,trial));
    hold on
    plot(a.(block_name{q(i)}).traj(:,7,trial),a.(block_name{q(i)}).traj(:,8,trial));
    pbaspect([1 1 1])
    axis([0.66 0.94 0.16 0.44])
    title(block_name{q(i)})
end