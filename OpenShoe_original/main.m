%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%> @file Main.m
%>                                                                           
%> @brief This is the skeleton file for demonstrating smoothing for
%> ZUPT-aided INS.
%> 
%> @details This is the skeleton file for demonstrating smoothing for 
%> ZUPT-aided INS. The data is first processed by a customary complementary
%> Kalman type of filter. Thereafter it is processed by a 3-pass smoothing
%> algorithms. The results are displayed overlaid on eachother in plots.
%> 
%>
%> @authors Isaac Skog, John-Olof Nilsson, David Simón Colomar
%> @copyright Copyright (c) 2011 OpenShoe, ISC License (open source)
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Loads the algorithm settings and the IMU data
disp('Loads the algorithm settings and the IMU data')
u=settings();

%% Run the zero-velocity detector 
disp('Runs the zero velocity detector')
[zupt T]=zero_velocity_detector(u);

%% Run the costumary Kalman filter
disp('Runs the customary Kalman filter')
[x_closed cov_closed] = ZUPTaidedINS(u, zupt);

%% Run the smoothing filter
disp('Runs the smoothing filter')
[x_smoothed, cov_smoothed, segment, P, P_smooth, dx, dx_smooth] = smoothed_ZUPTaidedINS(u, zupt);

%% View the result 
disp('Views the data')
view_data_smoothed;


