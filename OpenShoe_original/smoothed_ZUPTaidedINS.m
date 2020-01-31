%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%> @file smoothed_ZUPTaidedINS.m
%>
%> @brief This file contains all the functions needed to implement an 
%> open-loop zero-velocity aided inertial navigation system and smooth the 
%> estimated trajectory
%>
%> @details This file contains all the functions needed to implement an
%> open-loop zero-velocity aided inertial navigation system, given a set
%> of IMU data, initial closed loop estimations of the state and 
%> correlation matrix and a vector that indicates when the system has
%> zero-velocity.
%>
%> @authors David Simon Colomar, John-Olof Nilsson, Isaac Skog
%> @copyright Copyright (c) 2011 OpenShoe, ISC License (open source)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% MAINFUNCTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function [x_smoothed cov_smoothed] = smoothed_ZUPTaidedINS(u, zupt, x_closed, P_closed, quat_closed)
%
%> @brief Function that runs the open-loop zero-velocity aided INS Kalman
%> filter algorithm, applying an RTS smoothing filter
%>
%> @details Function that runs the open-loop zero-velocity aided INS Kalman
%> filter algorithm, applying an RTS  smoothing filter. All settings for
%> the filter are done in setting.m. 
%>
%> @param[out]  x               Matrix with the estimated navigation states.        
%> @param[out]  cov             Matrix with the diagonal elements of the state covariance matrices.       
%> @param[in]   u               The IMU data vector.     
%> @param[in]   zupt            Vector with the decisions of the zero-velocity.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x cov_smooth seg P P_smooth dx dx_smooth] = smoothed_ZUPTaidedINS(u, zupt)

global simdata;

% Length of data
N = length(zupt);

% Allocate state vectors
x = zeros(9, N);
quat = zeros(4,N);
dx = zeros(9, N);
dx_timeupd = zeros(9, N);
dx_smooth = zeros(9, N);

% Allocate covariance matrices
cov = zeros(9, N);
cov_smooth = zeros(9,N);
P = zeros(9,9,N);
P_timeupd = zeros(9,9,N);
P_smooth = zeros(9,9,N);
F = zeros(9,9,N);

% Constant matrices
[Q R H] = init_filter;
Id = eye(9);

% Initialize covariance matrix
P(1:3,1:3,1)=diag(simdata.sigma_initial_pos.^2);
P(4:6,4:6,1)=diag(simdata.sigma_initial_vel.^2);
P(7:9,7:9,1)=diag(simdata.sigma_initial_att.^2);
cov(:, 1) = diag(P(:,:,1));

% Initialize state vectors
[x(:, 1) quat(:,1)] = init_Nav_eq(u);

% Segment start/stop place holders
seg_start = 2;
seg_end = N;
seg = 1;
% Segment counter
c = 0;
while (1)

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Forward Kalman filter %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    for n = seg_start:seg_end

        
        %%%%%%%%%%%%%%%
        % Time update %
        %%%%%%%%%%%%%%%
        
        % Update mechanization equations
        [x(:, n) quat(:, n)] = Navigation_equations(x(:,n-1), u(:,n), quat(:,n-1));
        
        % Calculate state matrices
        [F(:, :, n) G] = state_matrix(quat(:, n), u(:,n));
        
        % Propagate errors
        dx(:,n) = F(:, :, n) * dx(:,n-1);
        
        % Update covariance matrices
        P(:, :, n) = F(:, :, n)*P(:, :, n - 1)*F(:, :, n)' + G*Q*G';
        
        % Save updated values for the smoothing
        dx_timeupd(:, n) = dx(:, n);
        P_timeupd(:, :, n) = P(:, :, n);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %      Zero-velocity update      %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if (zupt(n) == true)
            
            % Calculate the Kalman filter gain
            K=(P(:, :, n)*H')/(H*P(:, :, n)*H'+R);
            
            % Update deviation estimate
            dx(:, n)= dx(:, n) - K*(dx(4:6, n) - x(4:6,n));
            
            % Update the filter state covariance matrix P.
            P(:, :, n) = (Id - K*H)*P(:, :, n);
            
        end
        
        % Make sure the filter state covariance matrix is symmetric.
        P(:, :, n) = (P(:, :, n) + P(:, :, n)')/2;
        
        % Store the diagonal of the state covariance matrix P.
        cov(:, n) = diag(P(:, :, n));
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Segmentation decision %
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if c > 0
            c = c + 1;
        end
        if (sum(cov(4:6, n - 1)) > 0.1e-3) && (sum(cov(4:6, n)) < 0.1e-3) && (c == 0)
            c = 1;
        end
        if c == 30
            seg_end = n;
            c = 0;
            break;
        end
        
    end

    
    %%%%%%%%%%%%%%%%%
    % RTS smoothing %
    %%%%%%%%%%%%%%%%%
    
    % Initialize smoothing variables
    dx_smooth(:,seg_end)=dx(:,seg_end);
    P_smooth(:,:,seg_end)=P(:,:,seg_end);
    cov_smooth(:,seg_end) = diag(P_smooth(:,:,seg_end));
    
    for n = (seg_end-1):-1:seg_start
        
        % Variables needed for update equations
        A = P(:, :, n)*F(:, :, n)'/P_timeupd(:, :, n + 1);
        
        % Update state
        dx_smooth(:, n) = dx(:, n) + A*(dx_smooth(:, n + 1) - dx_timeupd(:, n + 1));
        
        % Update covariance
        P_smooth(:, :, n) = P(:, :, n) + A*(P_smooth(:, :, n+1) - P_timeupd(:, :, n + 1))*A';
        
        % Make sure the filter state covariance matrix is symmetric.
        P_smooth(:, :, n) = (P_smooth(:, :, n) + P_smooth(:, :, n)')/2;
        
        % Store the diagonal of the state covariance matrix P.
        cov_smooth(:,n) = diag(P_smooth(:, :, n));
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Internal state compensation %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for n = seg_start:seg_end

        [x(:,n) quat(:, n)]=comp_internal_states(x(:,n),- dx_smooth(:, n),quat(:, n));
    
    end
    
    
    %%%%%%%%%%%%%%%%%
    % Miscellaneous %
    %%%%%%%%%%%%%%%%%
    
    % Zero out last dx-value which will be used in next iteration
    dx(:,seg_end)=zeros(1,9);
    
    % Zero out these values to avoid linearization problems
    P(1:2,9,seg_end)=0;
    P(9,1:2,seg_end)=0;

    % Save segmentation points
    seg = [seg seg_end];
    
    % Check if end of data or correct segment place holders
    if seg_end~=N
        seg_start = seg_end + 1;
        seg_end = N;
    else
        % End of data -> break loop
        break
    end

    % Go back and process new segment
end

end


%% AUXILIAR FUNCTIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  funtion [x quat]=init_Nav_eq(u)
%
%> @brief Function that calculates the initial state of the navigation
%> equations.
%>
%> @details Function that calculates the initial state of the navigation
%> equations. That is, it does a simple initial alignment of the navigation
%> system, where the roll and pitch of the system is estimated from the
%> 20 first accelerometer readings. All other states are set according to
%> the information given in the function "settings.m".
%>
%> @param[out]  x     Initial navigation state vector.
%> @param[out]  quat  Quaternion vector, representating the initial attitude of the platform.
%> @param[in]   u     Matrix with the IMU data.
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x quat]=init_Nav_eq(u)

global simdata;

% Under the assumption that the system is stationary during the first 20
% samples, the initial roll and pitch is calculate from the 20 first
% accelerometer readings.
f_u=mean(u(1,1:20));
f_v=mean(u(2,1:20));
f_w=mean(u(3,1:20));

roll=atan2(-f_v,-f_w);
pitch=atan2(f_u,sqrt(f_v^2+f_w^2));


% Set the attitude vector
attitude=[roll pitch simdata.init_heading]';

% Calculate quaternion corresponing to the initial attitude
Rb2t=Rt2b(attitude)';
quat=dcm2q(Rb2t);

% Set the initial state vector
x=zeros(9,1);
x(1:3,1)=simdata.init_pos;
x(7:9,1)=attitude;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  funtion [Q R H]=init_filter()
%
%> @brief Function that initializes the Kalman filter.
%>
%> @details Function that initializes the Kalman filter. That is, the
%> function generates the initial covariance matrix P, the process noise
%> covariance matrix Q, the measurement noise covariance matrix R, and
%> observation matrix H, based upon the settings defined in the function
%> settings.m
%>
%> @param[out]   Q     Process noise covariance matrix.
%> @param[out]   R     Measurement noise covariance matrix.
%> @param[out]   H     Measurement observation matrix.
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Q R H]=init_filter

global simdata;

% Only the standard errors included
    
% Process noise covariance matrix
Q=zeros(6);
    
% Observation matrix
H=zeros(3,9);

% General values for the observation matrix H
H(1:3,4:6)=eye(3);

% General values for the process noise covariance matrix Q
Q(1:3,1:3)=diag(simdata.sigma_acc.^2);
Q(4:6,4:6)=diag(simdata.sigma_gyro.^2);

% General values for the measurement noise matrix R
R=diag(simdata.sigma_vel.^2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function q=dcm2q(R)
%
%>
%> @brief Function that converts a directional cosine matrix (rotation 
%> matrix) in to a quaternion vector. 
%>
%> @param[out]    q      Quaternion vector.
%> @param[in]     R      Rotation matrix.
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q=dcm2q(R)

T = 1 + R(1,1) + R(2,2) + R(3,3);

if T > 10^-8
    
    S = 0.5 / sqrt(T);
    qw = 0.25 / S;
    qx = ( R(3,2) - R(2,3) ) * S;
    qy = ( R(1,3) - R(3,1) ) * S;
    qz = ( R(2,1) - R(1,2) ) * S;

else
    
    if (R(1,1) > R(2,2)) && (R(1,1) > R(3,3))
        
        S = sqrt( 1 + R(1,1) - R(2,2) - R(3,3)) * 2; % S=4*qx
        qw = (R(3,2) - R(2,3)) / S;
        qx = 0.25 * S;
        qy = (R(1,2) + R(2,1)) / S;
        qz = (R(1,3) + R(3,1)) / S;
        
    elseif (R(2,2) > R(3,3))
        
        S = sqrt( 1 + R(2,2) - R(1,1) - R(3,3) ) * 2; %S=4*qy
        qw = (R(1,3) - R(3,1)) / S;
        qx = (R(1,2) + R(2,1)) / S;
        qy = 0.25 * S;
        qz = (R(2,3) + R(3,2)) / S;

    else
        
        S = sqrt( 1 + R(3,3) - R(1,1) - R(2,2) ) * 2; % S=4*qz
        qw = (R(2,1) - R(1,2)) / S;
        qx = (R(1,3) + R(3,1)) / S;
        qy = (R(2,3) + R(3,2)) / S;
        qz = 0.25 * S;

    end

end

%Store in vector
q = [qx qy qz qw]';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function R=q2dcm(q)
%
%>
%> @brief Function that converts a  quaternion vector to a directional
%> cosine matrix (rotation matrix) 
%>
%> @param[out]   R      Rotation matrix.
%> @param[in]    q      Quaternion vector.
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R=q2dcm(q)

p=zeros(6,1);

p(1:4)=q.^2;

p(5)=p(2)+p(3);

if p(1)+p(4)+p(5)~=0
   p(6)=2/(p(1)+p(4)+p(5)); 
else
   p(6)=0;
end


R(1,1)=1-p(6)*p(5);
R(2,2)=1-p(6)*(p(1)+p(3));
R(3,3)=1-p(6)*(p(1)+p(2));

p(1)=p(6)*q(1); 
p(2)=p(6)*q(2);
p(5)=p(6)*q(3)*q(4);
p(6)=p(1)*q(2);

R(1,2)=p(6)-p(5);
R(2,1)=p(6)+p(5);

p(5)=p(2)*q(4);
p(6)=p(1)*q(3);

R(1,3)=p(6)+p(5);
R(3,1)=p(6)-p(5);

p(5)=p(1)*q(4);
p(6)=p(2)*q(3);

R(2,3)=p(6)-p(5);
R(3,2)=p(6)+p(5);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function R=Rt2b(ang)
%
%>
%> @brief Function that calculates the rotation matrix for rotating a 
%> vector from coordinate frame t to the coordinate frame b, given a
%> vector of Euler angles.
%>
%> @param[out]  R      Rotation matrix.
%> @param[in]   ang    Euler angles [roll,pitch,heading]
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R=Rt2b(ang)


cr=cos(ang(1));
sr=sin(ang(1));

cp=cos(ang(2));
sp=sin(ang(2));

cy=cos(ang(3));
sy=sin(ang(3));

R=[cy*cp sy*cp -sp; 
    -sy*cr+cy*sp*sr cy*cr+sy*sp*sr cp*sr; 
    sy*sr+cy*sp*cr -cy*sr+sy*sp*cr cp*cr];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  funtion [y,q]=Navigation_equations(x,u,q)
%
%> @brief The mechanized navigation equations of the inertial navigation 
%> system. 
%>
%> @details The mechanized navigation equations of the inertial navigation 
%> system. That is, the function takes the old state (position, velocity, 
%> and attitude) of the navigation system, togheter with the current IMU
%> data measurements (specific force, angular rates), and calculates 
%> the current state of the navigation system.  
%>
%> @Note The mechanization of the navigation equations that has been 
%> implemented is very simple, and several higher order terms has been 
%> neglected. Therefore, this mechanization of the navigation equations 
%> should only be used in systems using low-cost sensor and where only 
%> moderate velocities can be expected. 
%>
%> @param[out]   y     New navigation state [position,velocity, attitude (euler angles].
%> @param[out]   q     New quaternions
%> @param[in]    x     Old navigation state
%> @param[in]    u     IMU data [specific force, angular rates].
%> @param[in]    q     Old quaternions
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y,q]=Navigation_equations(x,u,q)

global simdata;

% Allocate memmory for the output vector
y=zeros(size(x));

% Get sampling period of the system
Ts=simdata.Ts;

%*************************************************************************%
% Update the quaternion vector "q"  given the angular rate measurements.
%*************************************************************************%

w_tb=u(4:6);

P=w_tb(1)*Ts;
Q=w_tb(2)*Ts;
R=w_tb(3)*Ts;

OMEGA=zeros(4);
OMEGA(1,1:4)=0.5*[0 R -Q P];
OMEGA(2,1:4)=0.5*[-R 0 P Q];
OMEGA(3,1:4)=0.5*[Q -P 0 R];
OMEGA(4,1:4)=0.5*[-P -Q -R 0];

v=norm(w_tb)*Ts;

if v~=0
    q=(cos(v/2)*eye(4)+2/v*sin(v/2)*OMEGA )*q;
    q=q./norm(q);
end

%*************************************************************************%
% Use the update quaternion to get attitude of the navigation system in
% terms of Euler angles.
%*************************************************************************%

% Get the roll, pitch and yaw
Rb2t=q2dcm(q);
% roll
y(7)=atan2(Rb2t(3,2),Rb2t(3,3));

% pitch
y(8)= atan(- Rb2t(3,1)/sqrt(Rb2t(3, 2)^2 + Rb2t(3,3)^2));

%yaw
y(9)=atan2(Rb2t(2,1),Rb2t(1,1));


%*************************************************************************%
% Update position and velocity states using the measured specific force,
% and the newly calculated attitude.
%*************************************************************************%

% Gravity vector
g_t=[0 0 simdata.g]';

% Transform the specificforce vector into navigation coordinate frame.
f_t=q2dcm(q)*u(1:3);

% Subtract (add) the gravity, to obtain accelerations in navigation
% coordinat system.
acc_t=f_t+g_t;

% State space model matrices
A=eye(6);
A(1,4)=Ts;
A(2,5)=Ts;
A(3,6)=Ts;

B=[(Ts^2)/2*eye(3);Ts*eye(3)];

% Update the position and velocity estimates.
y(1:6)=A*x(1:6)+B*acc_t;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  funtion [F G]=state_matrix(q,u)
%
%> @brief Function for calculating the state transition matrix F and 
%> the process noise gain matrix G. 
%>
%> @details Function for calculating the state transition matrix F and 
%> the process noise gain matrix G, given the current orientation of 
%> the platform and the specific force vector.  
%>
%> @param[out]   F     State transition matrix.
%> @param[out]   G     Process noise gain matrix.
%> @param[in]    u     IMU data [specific force, angular rates].
%> @param[in]    q     Old quaternions
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F G]=state_matrix(q,u)

global simdata

% Convert quaternion to a rotation matrix
Rb2t=q2dcm(q);

% Transform measured force to force in
% the navigation coordinate system.
f_t=Rb2t*u(1:3);

% Create a ske symmetric matrix of the specific fore vector
% This is a vectorial product matrix
St=[0 -f_t(3) f_t(2); f_t(3) 0 -f_t(1); -f_t(2) f_t(1) 0];

% Zero matrix
O=zeros(3);

% Identity matrix
I=eye(3);

% Transition matrix
    Fc=[O I O;
        O O St;
        O O O];
    
% Noise gain matrix
    Gc=[O O; Rb2t O; O -Rb2t];



% Approximation of the discret time state transition matrices
F=eye(size(Fc))+simdata.Ts*Fc;
G=simdata.Ts*Gc;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function x_out = comp_internal_states(x_in, dx, q_in)
%
%> @brief Function that corrects the estimated navigation states with 
%> the by the Kalman filter estimated  system perturbations (errors) and
%> the open loop estimation of the signal
%>
%> @details Function that corrects the estimated navigation states with 
%> the by the Kalman filter estimated system perturbations (errors). That
%> is, the current position an velocity estimates of the navigation 
%> platform is corrected by adding the estimated system perturbations to 
%> these states. To correct the orientation state (Euler angles), the
%> quaternion vector are first converted into a rotation matrix, which then
%> is corrected using the estimated orientation perturbations. The
%> corrected rotation matrix is then transformed back into a quaternion
%> vector, as well as the equivalent vector of Euler angles.         
%>
%> @param[out]   x_out     Corrected (posteriori) navigation state vector.
%> @param[in]    x_in      A priori estimated navigation state vector.
%> @param[in]    dx        Vector of system perturbations calculated by Kalman filter
%> @param[in]    q_in      A priori estimated quaternion vector.
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_out q_out]=comp_internal_states(x_in,dx,q_in)

% Convert quaternion to a rotation matrix
R=q2dcm(q_in);

% Correct the state vector
x_out=x_in+dx;

% Correct the rotation matrics
epsilon=dx(7:9);
OMEGA=[0 -epsilon(3) epsilon(2); epsilon(3) 0 -epsilon(1); -epsilon(2) epsilon(1) 0];
R=(eye(3)-OMEGA)*R;


% Get the corrected roll, pitch and heading from the corrected rotation
% matrix
x_out(7)=atan2(R(3,2),R(3,3));
x_out(8)=atan(-R(3,1)/sqrt(R(3,2)^2+ R(3,3)^2));
x_out(9)=atan2(R(2,1),R(1,1));

% Calculte the corrected quaternions
q_out=dcm2q(R);

end
