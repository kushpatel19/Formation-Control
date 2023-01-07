% First Control Law
%%
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
n=4;            % No of Robots
xc = 5*sin(t);       % X coordinate of virtual frame
yc = 5*cos(t);       % Y coordinate of virtual frame
theta_c = 0;    % Orientation of virtual frame w.r.t inertial frame
zeta = [xc yc theta_c];  % State of virtual coordinate frame
G = [[1 2 0 2];          % Interaction Topology (represent that how other
     [2 1 2 0];          % robots are connected to each others)
     [0 2 1 2];
     [2 0 2 1]];
g = size(G);    % Size of interation topology matrix
% disp(g)
alpha = 1;      % proportional gain (positive scalar)
gama = 2;       % proportional gain (positive scalar)
v = diff(xc)/dt;    % Linear velocity of robot
% disp(length(v))
% disp(length(t))
omega_c = diff(theta_c);   % Angular velocity of robot
rf_d = [[-1 1];            % Desired deviation w.r.t virtual coordinate frame
        [1 1];
        [1 -1];
        [-1 -1]];
% r_d = rf_d;    
r0=  [[-3 2];            % Initial position of all robots
    [0 0];
    [7 9];
    [-3 8]];
R_d = zeros(2*(length(t)-1), g(1));   % Desired position-time data storage
Xd = zeros(g(1),(tf - t0)/(dt));      % Desired X-position data storage
Yd = zeros(size(Xd));                 % Desired Y-position data storage
X = zeros(2*(length(t)-1)+2, g(1));   % Real position-time data storage
X(1:2,:) = r0';                        % Initialize the real position at t=t0

%% State Transition Matrix 'R_d'

R = [cos(theta_c) -sin(theta_c);sin(theta_c) cos(theta_c)]; % Rotation Matrix
for i= 1:(length(t)-1)
%   t = t(i);
%   vc = v(i);
    cf = [xc(i);yc(i)];   % State of virtual frame
    for j=1:g(1)
        C = [rf_d(j,1);rf_d(j,2)];
%       disp(C)
        D = cf + (R*C);   
%       disp(D)
        R_d(2*i-1:2*i,j) = D;
        Xd(j,i) = D(1,1);
        Yd(j,i) = D(2,1);
%       disp(j)
%       disp(Xd(j,i))
    end
end
% disp(Xd)
% disp(Yd)
% disp(size(Xd))
% disp(size(Yd))
Xd_dot = (diff(Xd'))'/dt;
Xd_dot = [Xd_dot Xd_dot(:,end)];   % Desired X-velocity data storage
Yd_dot = (diff(Yd'))'/dt;
Yd_dot = [Yd_dot Yd_dot(:,end)];   % Desired Y-velocity data storage
% disp(size(Xd_dot))
% disp(size(Yd_dot))
% disp(Xd_dot)
% disp(Yd_dot)

%% Consensus Algorithm

for i = 1:(length(t)-1)
    for j=1:g(1)
        p = [Xd_dot(j,i); Yd_dot(j,i)];
        q = [X(2*i-1,j); X(2*i,j)];
        r = [R_d(2*i-1,j); R_d(2*i,j)];
        u = p - (alpha*(q-r));        % Control Law (Kinematic Equation)
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
