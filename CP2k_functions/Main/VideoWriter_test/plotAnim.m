clc;
clear all;

v= VideoWriter('exp2dvideo.avi');
open(v);


load('torque');
for i=1:2:40
    for j = 1:5000
        tqAnim(i,j)=imTorque(1,j+(i-1)*2500);
        tqAnim(i+1,j)=imTorque(2,j+(i-1)*2500);
    end
end


figure;
hold on;
for i=1:20
    plot(tqAnim((i-1)*2+1,:),tqAnim((i-1)*2+2,:));
    axis([0 1 -40 40]);
    xlabel('Time in sec')
    ylabel('Torque in Nm')
    title('Plot of developed Torque')
    set(findobj(gca,'Type','line','Color',[0 0 1]),...
    'Color','blue',...
    'LineWidth',1.2)
    F(i) = getframe(gcf);
end
hold off;

close(v);
%movie2avi(F,'torque','compression','None')
