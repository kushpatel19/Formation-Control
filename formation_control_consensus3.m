%%
% Third Control Law
clc;
clear all;
close all;
clear;

%% Info
% Consensus Algorithm (Distributed Formatation Control Architecture)
% written by Kush Patel        [Copyright :)]

%% Initial Conditions & Parameters
t0= 0;          % Initial time
tf = 10;        % Final time
dt = 0.1;       % Time step
N = (tf-t0)/dt; % Iterations
t = t0:dt:tf;
n=4;            % No. of Robots
xc = (t);       % X coordinate of virtual frame
yc = 0.1*(t.*t);       % Y coordinate of virtual frame
xc_0 = xc(1);   % X coordinate of virtual frame at t=t0
yc_0 = xc(1);   % Y coordinate of virtual frame at t=t0
theta_c = 0;    % Orientation of virtual frame w.r.t inertial frame
% zeta = [xc yc theta_c];  % State of virtual coordinate frame
zeta_0 = [xc_0, yc_0];     % State of virtual coordinate frame at t=t0
G = [[1 2 0 2 2];          % Interaction Topology (represent that how other
     [2 1 2 0 2];          % robots are connected to each others)
     [0 2 1 2 2];
     [2 0 2 1 2];          % 1 extra row and column due to virtual robot
     [2 2 2 2 1]];
g = size(G);    % Size of interation topology matrix  (n+1,n+1)
alpha = 1;      % proportional gain (positive scalar)
gama = 2;       % proportional gain (positive scalar)
v = diff(xc)/dt;    % Linear velocity of robot
omega_c = diff(theta_c);   % Angular velocity of robot
rf_d = [[-1 1];            % Desired position of the formation 
        [1 1];             % w.r.t virtual coordinate frame
        [1 -1];
        [-1 -1]];
%  rf_d = [[0 -1];            % Desired position w.r.t virtual coordinate frame
%          [0 -2];
%          [0 -3];
%          [0 -4]];
% r_d = rf_d;    
r0=  [[-4 2];            % Initial position of all robots
      [2 3];
      [2 -4];
      [-4 -2]];

R_d = zeros(2*(length(t)-1), n);    % Desired position-time data storage
Xd = zeros(n,N);                    % Desired X-position data storage
Yd = zeros(size(Xd));               % Desired Y-position data storage
X = zeros(2*(length(t)-1)+2, n);    % Real position-time data storage
X(1:2,:) = r0';                     % Initialize the real position at t=t0
ZETA = zeros(2*(length(t)-1)+2, n); % Real position-time data storage of virtial center
ZETA(1:2,:) = r0';                  % Initialize the real position at t=t0
zeta_d = zeros(2,(length(t)));      % Desired position-time data storage of virtual center
zeta_d(:,1) = zeta_0';              % Initialize the real position at t=t0
zeta = zeros(2*(length(t)),n);      % Real position-time data storage of virtial center
zeta(1,:) = sum(R_d(1,:))/4;        % Initialize the X-coordinate of real position at t=t0
zeta(2,:) = sum(R_d(2,:))/4;        % Initialize the Y-coordinate of real position at t=t0

%% State Transition Matrix 'R_d' & Desired position-time Matrix of Virtual Robot/Center

R = [cos(theta_c) -sin(theta_c);sin(theta_c) cos(theta_c)]; % Rotation Matrix
for i= 1:(length(t)-1)
%     t = t(i);
%     vc = v(i);
    cf = [xc(i);yc(i)];   % State of virtual frame
    for j=1:n
        C = [rf_d(j,1);rf_d(j,2)];
        D = cf + (R*C);   
        R_d(2*i-1:2*i,j) = D;
        Xd(j,i) = D(1,1);
        Yd(j,i) = D(2,1);
    end
end

for i = 1:(length(t)-1)
    zeta_d(1,i+1) = sum(R_d(2*i-1,:))/4;  % X-coordinate of desired zeta
    zeta_d(2,i+1) = sum(R_d(2*i,:))/4;    % Y-coordinate of desired zeta 
end
zeta_d_dot = (diff(zeta_d'))'/dt;         % derivative of desired zeta
zeta_d_dot = [zeta_d_dot zeta_d_dot(:,end)];

%% Consensus Algorithm

% Consensus in zeta
for i = 1:(length(t)-1)
    for j = 1:n
        [aa,bb] = sensor(ZETA(2*i-1,j),ZETA(2*i,j));
        [cc,dd] = sensor(ZETA(2*i-1,j),ZETA(2*i,j));
        [ee,ff] = sensor(ZETA(2*i-1,j),ZETA(2*i,j));
        [gg,hh] = sensor(ZETA(2*i-1,j),ZETA(2*i,j));
        zeta(2*i+1,j) = (aa+cc+ee+gg)/4;
        zeta(2*i+2,j) = (bb+dd+ff+hh)/4;
    end
    zeta_dot = (zeta(2*i+1:2*i+2,:) - zeta(2*i-1:2*i,:))/dt;
    for m = 1:n
        aa = zeta_d_dot(:,i+1);
        bb = [zeta(2*i+1,m); zeta(2*i+2,m)];
        cc = zeta_d(:,i+1);
        zz = 0;
        for p = 1:n
            dd = zeta_dot(:,p);
            ee = [zeta(2*i+1,p); zeta(2*i+2,p)];
            zz = zz + G(m,p)*(dd - gama*(bb - ee));
        end
        ni = sum(G(m,:));
        u = (zz/ni) + (1/ni)*G(m,n+1)*(aa - gama*(bb - cc)); % Control Law for zeta
        
        ZETAdot = u;
        ZETA(2*i+1:2*i+2,m) = ZETA(2*i-1:2*i,m) + ZETAdot*dt; % Euler Integration
    end
end
% Final (average) zeta
zzz = zeros(2,length(t));
for i = 1:length(t)
    zzz(1,i) = sum(ZETA(2*i-1,:))/4;
    zzz(2,i) = sum(ZETA(2*i,:))/4;
end
% Updating desired valued based on new zeta
for i= 1:(length(t)-1)
%   t = t(i);
%   vc = v(i);
    cf = [zzz(1,i);zzz(2,i)];   % State of virtual robot/center
    for j=1:n
        C = [rf_d(j,1);rf_d(j,2)];
        D = cf + (R*C);   
        R_d(2*i-1:2*i,j) = D;
        Xd(j,i) = D(1,1);
        Yd(j,i) = D(2,1);
    end
end

Xd_dot = (diff(Xd'))'/dt;
Xd_dot = [Xd_dot Xd_dot(:,end)];   % Desired X-velocity data storage
Yd_dot = (diff(Yd'))'/dt;
Yd_dot = [Yd_dot Yd_dot(:,end)];   % Desired Y-velocity data storage

% Consesus in position
for i = 1:(length(t)-1)
    for j=1:n
        p = [Xd_dot(j,i); Yd_dot(j,i)];
        q = [X(2*i-1,j); X(2*i,j)];
        r = [R_d(2*i-1,j); R_d(2*i,j)];
        l=0;
        for k = 1:n
            z = [X(2*i-1,k); X(2*i,k)];
            w = [R_d(2*i-1,k); R_d(2*i,k)];
            l = l + G(j,k)*((q -r) - (z -w));
        end
        u = p - (alpha*(q-r)) - l;        % Control Law (Kinematic Equation)
        Xdot = u;
        X(2*i+1:2*i+2,j) = X(2*i-1:2*i,j) + Xdot*dt;   % Euler Integration
    end
end

%% plot without animation

for i=1:2:2*(length(t)-1)
    plot([X(i,:),X(i,1)],[X(i+1,:),X(i+1,1)],'-')
    xlim([min(min(X))-1,max(max(X))+1]);
    ylim([min(min(X))-1,max(max(X))+1]);
    hold on
end

%% Animate and save video
figure()

%%% create list of colors %%%
C = {'k','b','r','g'}; % list of 4 colors
cl= 1;
%% plot initial positions %%%

% for i = 1:2:2*length(t)-1
%         hold on
%         plot([X(i,:),X(i,1)],[X(i+1,:),X(i+1,1)],'LineWidth',2,'color',C{cl},'marker','o');
%         cl = cl + 1;
%         if cl == 5
%             cl = 1;
%         end
% end
c =1;
%%% create a video writer object %%%
myVideo = VideoWriter('myVideoFile','MPEG-4'); %open video file
myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
open(myVideo)
u = uicontrol('Style','slider',...% create slider indicating time remaining
    'Min',1,'Max',201,'Value',1);

%%% animate using drawnow %%%
    for i = 1:2:2*length(t)-1
        plot([X(i,:),X(i,1)],[X(i+1,:),X(i+1,1)],'LineWidth',1,'color',C{1},'marker','o');
        xlim([min(min(X))-1,max(max(X))+1]);
        ylim([min(min(X))-1,max(max(X))+1]);
        drawnow limitrate;
%       pause(0.01)
        frame = getframe(gcf);
        writeVideo(myVideo, frame);    % save video as MP4
        clf;
        c = c + 1;
        if c == 4
            c = 1;
        end
    end
%     u.Value = i; % update slider
close(myVideo)

%% Function of sensor
function [a,b]=sensor(x,y)
    a = normrnd(x,0.1);
    b = normrnd(y,0.1);
end