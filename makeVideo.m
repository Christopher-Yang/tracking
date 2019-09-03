
fhandle = figure(1); clf; hold on
    set(fhandle, 'Position', [50, 50, 800, 800]); % set size and loction on screen
    set(fhandle, 'Color','w') % set background color to white
    
v = VideoWriter('test.avi','Motion JPEG AVI');
v.FrameRate = 130;
open(v);
figure(1); clf; hold on

for t=1:5201
    cla
    axis off
    plot(dat.target.x_pos(1:t),dat.target.y_pos(1:t),'k')
    plot(dat.Rhand.x_pos(1:t),dat.Rhand.y_pos(1:t),'r')
    plot(dat.target.x_pos(t),dat.target.y_pos(t),'ko','markersize',12,'markerfacecolor',[1 1 0],'linewidth',2)
    plot(dat.Rhand.x_pos(t),dat.Rhand.y_pos(t),'ko','markersize',8,'markerfacecolor',[1 0 0],'linewidth',1.5)
    
    plot(t_onFreq(1:t,1)+0.25,t_onFreq(1:t,2),'k')
    plot(h_onFreq(1:t,1)+0.25,h_onFreq(1:t,2),'r')
    plot(t_onFreq(t,1)+0.25,t_onFreq(t,2),'ko','markersize',12,'markerfacecolor',[1 1 0],'linewidth',2)
    plot(h_onFreq(t,1)+0.25,h_onFreq(t,2),'ko','markersize',8,'markerfacecolor',[1 0 0],'linewidth',1.5)

    plot(t_onFreq(1:t,1) - h_onFreq(1:t,1)+0.5,t_offFreq(1:t,2) - h_onFreq(1:t,2),'k')
    plot(h_offFreq(1:t,1)+0.5,h_offFreq(1:t,2),'r')
    plot(t_onFreq(t,1) - h_onFreq(t,1)+0.5,t_offFreq(t,2) - h_onFreq(t,2),'ko','markersize',12,'markerfacecolor',[1 1 0],'linewidth',2)
    plot(h_offFreq(t,1)+0.5,h_offFreq(t,2),'ko','markersize',8,'markerfacecolor',[1 0 0],'linewidth',1.5)
    axis equal
    axis([-.12 .62 -.12 .12])
    
    frame = getframe;
    writeVideo(v,frame);
end
close(v);