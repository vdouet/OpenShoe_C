%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%> @file view_data_smoothed.m  
%>
%> @brief Script for plotting the data from the zero-velocity aided inertial
%> navigations system.
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global simdata;


% Close all windows
close all

% Generate a vector with the time scale
N=length(cov_closed);
t=0:simdata.Ts:(N-1)*simdata.Ts;

%% Plot the IMU data
figure(1)
clf
subplot(2,1,1)
plot(t,u(1:3,:)')
xlabel('time [s]')
ylabel('Specific force [m/s^2]')
title('Specific force (accelerometer) measurements')
legend('x-axis','y-axis','z-axis')
box on
grid on

subplot(2,1,2)
plot(t,u(4:6,:)'*180/pi)
xlabel('time [s]')
ylabel('Angular rate  [deg/s]')
title('Angular rate measurements')
legend('x-axis','y-axis','z-axis')
box on
grid on


%% Plot the trajectory in the horizontal plane
figure(2)
clf
plot(x_closed(2, :),x_closed(1, :))
hold;
plot(x_smoothed(2, :), x_smoothed(1, :), 'g')
plot(x_closed(2,1),x_closed(1,1),'rs')
title('Trajectory')
legend('Trajectory', 'Smoothed trajectory', 'Start point')
xlabel('x [m]')
ylabel('y [m]')
axis equal
grid on
box on

%% Plot the height profile, the speed and when ZUPTs were applied SMOOTHING

figure(3)
clf
subplot(3,1,1)
plot(t, -x_closed(3, :))
hold;
plot(t, -x_smoothed(3, :), 'g')
title('Heigth')
xlabel('time [s]')
ylabel('z [m]')
grid on
box on
legend('closed loop', 'smoothed')


subplot(3,1,2)
plot(t,sqrt(sum(x_closed(4:6, :).^2)))
hold;
plot(t, sqrt(sum(x_smoothed(4:6,:).^2)), 'g')
title('Speed')
xlabel('time [s]')
ylabel('|v| [m/s]')
grid on
box on
legend('closed loop', 'smoothed')


subplot(3,1,3)
stem(t,zupt)
title('Zupt applied')
xlabel('time [s]')
ylabel('on/off')
grid on
box on


%% Plot the attitude

figure(4)
clf
subplot(2, 1, 1)
plot(t,unwrap(x_closed(7:9, :)')*180/pi)
title('Attitude')
xlabel('time [s]')
ylabel('Angle [deg]')
legend('Roll','Pitch','Yaw')
grid on
box on


subplot(2, 1, 2)
plot(t, unwrap(x_smoothed(7:9, :)')*180/pi)
title('Attitude (smoothed)')
xlabel('time [s]')
ylabel('Angle [deg]')
legend('Roll','Pitch','Yaw')
grid on
box on


%% Plot the diagonal elements of the filter covariance matrices as a
%% function of time

figure(5)
clf

subplot(3,2,1)
plot(t,sqrt(cov_closed(1:3, :))')
title('Position covariance')
ylabel('sqrt(cov) [m]')
xlabel('time [s]')
legend('x-axis', 'y-axis','z-axis')
grid on
box on


subplot(3,2,3)
plot(t,sqrt(cov_closed(4:6, :))')
title('Velocity covariance')
ylabel('sqrt(cov) [m/s]')
xlabel('time [s]')
legend('x-axis', 'y-axis','z-axis')
grid on
box on

subplot(3,2,5)
plot(t,sqrt(cov_closed(7:9, :))'*180/pi)
title('Heading covariance')
ylabel('sqrt(cov) [deg]')
xlabel('time [s]')
legend('Roll', 'Pitch','Yaw')
grid on
box on

% The same but for the smoothed
subplot(3,2,2)
plot(t,sqrt(cov_smoothed(1:3, :))')
title('Smoothed position covariance')
ylabel('sqrt(cov) [m]')
xlabel('time [s]')
legend('x-axis', 'y-axis','z-axis')
grid on
box on


subplot(3,2,4)
plot(t,sqrt(cov_smoothed(4:6, :))')
title('Smoothed velocity covariance')
ylabel('sqrt(cov) [m/s]')
xlabel('time [s]')
legend('x-axis', 'y-axis','z-axis')
grid on
box on

subplot(3,2,6)
plot(t,sqrt(cov_smoothed(7:9,:))'*180/pi)
title('Smoothed heading covariance')
ylabel('sqrt(cov) [deg]')
xlabel('time [s]')
legend('Roll', 'Pitch','Yaw')
grid on
box on

figure(6), hold on
plot(dx_smooth(1,:),'b')
plot(dx_smooth(2,:),'g')
plot(dx(1,:),'r')
plot(dx(2,:),'k')
title('Corrections xy')
ylabel('Position corrections [m]')
xlabel('Sample nr')
grid on
box on
legend('smoothed dx [m]','smoothed dy [m]','dx [m]','dy [m]')

figure(7), hold on
first_segment=5;
nr_segments_to_plot=15;
for n=first_segment:(nr_segments_to_plot+first_segment)
    % Cut out segment (smoothed)
    x_seg = x_smoothed(1:3,segment(n):segment(n+1));
    x_seg(1,:)=x_seg(1,:)-x_smoothed(1,segment(n));
    x_seg(2,:)=x_seg(2,:)-x_smoothed(2,segment(n));
    x_seg(3,:)=x_seg(3,:)-x_smoothed(3,segment(n));
    % Remove planar rotation over step
    R_seg = planerot(x_seg(1:2,end));
    x_seg(1:2,:)=R_seg*x_seg(1:2,:);
    % Remove difference in height
    x_seg(3,:) = x_seg(3,:) -linspace(x_seg(3,1),x_seg(3,end),size(x_seg,2));
    % Plot smoothed step
    plot3(x_seg(2,:)+n-first_segment,x_seg(1,:),-x_seg(3,:),'r');
    % Cut out segment (unsmoohted)
    x_seg = x_closed(1:3,segment(n):segment(n+1));
    x_seg(1,:)=x_seg(1,:)-x_closed(1,segment(n));
    x_seg(2,:)=x_seg(2,:)-x_closed(2,segment(n));
    x_seg(3,:)=x_seg(3,:)-x_closed(3,segment(n));
    % Remove planar rotation over step
    R_seg = planerot(x_seg(1:2,end));
    x_seg(1:2,:)=R_seg*x_seg(1:2,:);
    % Remove difference in height
    x_seg(3,:) = x_seg(3,:) -linspace(x_seg(3,1),x_seg(3,end),size(x_seg,2));
    % Plot unsmoothed step
    plot3(x_seg(2,:)+n-first_segment,x_seg(1,:),-x_seg(3,:),'b');
end
grid on
axis([-2 nr_segments_to_plot+2 -0.1 2 -0.02 0.25])
title('Steps'), xlabel('position x [m]'), ylabel('position x [m]'), zlabel('position x [m]')
set(gca,'CameraPosition',[74.7728   14.4268    1.6224])
set(gca,'CameraViewAngle',11.0958)
set(gca,'DataAspectRatio',[30 3.8 2])
zoom(1.4)
