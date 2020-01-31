
figure(7), hold on
for n=1:(length(segment)-1)
    x_seg = x_smoothed(1:3,segment(n):segment(n+1))-x_smoothed(1:3,segment(n));
    R_seg = planerot(x_seg(1:2,end));
    x_seg(1:2,:)=R_seg*x_seg(1:2,:);
    plot3(x_seg(1,:),x_seg(2,:),x_seg(3,:));
end