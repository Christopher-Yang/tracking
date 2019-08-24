subj = 3;
trial = 3;
a = d.(subj_rot{subj});
q = [1:2 5:6];

figure(1); clf
for i = 1:4
    subplot(1,4,i); hold on
    plot(a.(block_name{q(i)}).traj(1:650,1,trial),a.(block_name{q(i)}).traj(1:650,2,trial));
    plot(a.(block_name{q(i)}).traj(1:650,7,trial),a.(block_name{q(i)}).traj(1:650,8,trial));
    pbaspect([1 1 1])
    axis([0.66 0.94 0.16 0.44])
    title(block_name{q(i)})
end

subj = 7;
trial = 5;
a = d.(subj_rot_i{subj});
q = [1:2 5:6];

figure(2); clf
for i = 1:4
    subplot(1,4,i); hold on
    plot(a.(block_name{q(i)}).traj(1:650,1,trial),a.(block_name{q(i)}).traj(1:650,2,trial));
    plot(a.(block_name{q(i)}).traj(1:650,7,trial),a.(block_name{q(i)}).traj(1:650,8,trial));
    pbaspect([1 1 1])
    axis([0.66 0.94 0.16 0.44])
    title(block_name{q(i)})
end