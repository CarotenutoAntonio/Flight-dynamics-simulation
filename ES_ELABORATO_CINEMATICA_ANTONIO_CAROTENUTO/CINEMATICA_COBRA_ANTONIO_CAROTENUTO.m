clear all; close all; clc

%% Prepare to solve kinematics equations



% Initial conditions
psi0 = 0; theta0 = 0; phi0 = 0; % assi body orientato come assi Earth 
Z0 = 0;

%% Anonymous functions

% LOOPING
% three appropriate time histories of p, q, r
p_max = convangvel(0.0,'deg/s','rad/s'); % rad/s
q_max = convangvel(11,'deg/s','rad/s'); % rad/s 
r_max = convangvel(0,'deg/s','rad/s'); % rad/s

% tempo di simulazione
t_fin=28;

p = @(t) 0*t;

%velocità angolari
p = @(t) 0*t;

vBreakPointsQ(1,:) = [0 , 0.0179 , 0.0357 , 0.0714 , 0.1071 , 0.2857 ,...
                      0.3571 , 0.4643 , 0.7857 , 1.0714]*t_fin;

vBreakPointsQ(2,:) = [0 , 0.3907 , 0.9636 , 1.0417 , 2.1095 , 0.9063 ,...
                     -0.8292 , -1.4584 , -0.1042 , -0.0156]*q_max;
q = @(t) ...
    interp1( vBreakPointsQ(1,:), vBreakPointsQ(2,:), t, 'pchip' );

r = @(t) 0*t;

% three appropriate time histories of u, v, w
u0 = convvel(180 , 'km/h','m/s');
w0 = convvel(0 , 'km/h','m/s');

vBreakPointsU(1 ,:) = [0 , 0.0714 , 0.1786 , 0.3571 , 0.5 , 0.6429 ,...
                      0.8929 , 1.0714]*t_fin; 
vBreakPointsU(2 ,:) = [1 , 0.5 , 0.28 , 0.01 , 0.2 , 0.4 ,...
                      0.8 , 1.3]*u0;

u = @(t) ...
    interp1( vBreakPointsU(1,:), vBreakPointsU(2,:), t, 'pchip' );
% assume alpha = beta = 0
v = @(t) ...
    0;
w = @(t) ...
    0;

%% Anonymous functions
% RHS of quaternion components evolution equations, see eq. (2.27)
dQdt= @(t,Q) 0.5.*....
[ 0 -p(t) -q(t) -r(t) ;
p(t) 0 r(t) -q(t) ;
q(t) -r(t) 0 p(t) ;
r(t) q(t) -p(t) 0 ]*Q;

%% Solution of Gimbal equations
options = odeset( ...
    'RelTol', 1e-9, ...
    'AbsTol', 1e-9 ...
    );

Q0=angle2quat(psi0, theta0, phi0);

[vTime, vQ] = ode45(dQdt, [0 t_fin], Q0, options);


[Psi , Theta , Phi] = quat2angle(vQ);

vPhiThetaPsi=[Phi,Theta,Psi];


%% Integrate Navigation equations
% Time interpolation function for known Euler angles histories
fPhi = @(t) ...
    interp1(vTime,vPhiThetaPsi(:,1),t);
fTheta = @(t) ...
    interp1(vTime,vPhiThetaPsi(:,2),t);
fPsi = @(t) ...
    interp1(vTime,vPhiThetaPsi(:,3),t);
% RHS of navigation equations
dPosEdt = @(t,Pos) ...
    transpose(angle2dcm(fPsi(t),fTheta(t),fPhi(t),'ZYX'))*[u(t);v(t);w(t)]; %  + 0.*Pos;

%% Solution of navigation equations
options = odeset( ...
    'RelTol', 1e-3, ...
    'AbsTol', 1e-3*ones(3,1) ...
    );
PosE0 = [0;0;0];
[vTime2, vPosE] = ode45(dPosEdt, vTime, PosE0, options);
N = length(vPosE);
           
vXe = vPosE(:,1); vYe = vPosE(:,2); vZe = Z0 + vPosE(:,3);


%% Plots

% Euler angles time histories
figure(1)
subplot 121,
plot( ...
    vTime,convang(Psi,'rad','deg'),'-', ...
    vTime,convang(Theta,'rad','deg'),'--', ...
    vTime,convang(Phi,'rad','deg'),'-' ...
    )
legend('\phi','\theta','\psi')
xlabel('t (s)'); ylabel('(deg)')
title('Euler angles')
% set(gca,'fontname','cambria','fontsize',15)

% Angular velocity components
subplot 122
plot( ...
    vTime, convangvel(q(vTime),'rad/s','deg/s'), ...
    vBreakPointsQ(1,:), convangvel(vBreakPointsQ(2,:),'rad/s','deg/s'),'.');


    
% axis([0 20 -.3 .6])
legend('q(t)')
xlabel('t (s)'); ylabel('(deg/s)')
title('Angular velocity components in body axes')
% set(gca,'fontname','cambria','fontsize',15)


% CG coordinates time history
figure(2)
plot( ...
    vTime2,vXe, '-',...
    vTime2,vYe, '-.', ...
    vTime2,vZe, '--' ...
    )
% hold on, view(3)
legend('x_{G,E}(t)','y_{G,E}(t)','z_{G,E}(t)')
xlabel('t (s)'); ylabel('(m)');
title('CG coordinates in Earth axes');



%% Setup the figure/scene for 3D visualization
h_fig3 = figure(3);
grid on;
hold on;
light('Position',[1 0 -4],'Style','local');
% Trick to have Ze pointing downward and correct visualization
set(gca,'XDir','reverse');
set(gca,'ZDir','reverse');
daspect([1 1 1]);

%% Load aircraft shape
shapeScaleFactor = 25.0;
%shape = loadAircraftMAT('aircraft_pa24-250.mat', scale_factor);
shape = loadAircraftMAT('aircraft_mig29.mat', shapeScaleFactor);

mXYZe = [vPosE(:,1),vPosE(:,2),vPosE(:,3)+Z0];
mEulerAngles = [vPhiThetaPsi(:,3),vPhiThetaPsi(:,2),vPhiThetaPsi(:,1)];

%% Settings
% General settings
options.samples = [1,41,81,121,151:50:numel(vTime)]; %[1,40,80,120,160,200,250,300,320,numel(vTime)];
options.theView = [110 60 45];

% body axes settings
options.bodyAxes.show = true;
options.bodyAxes.magX = 1.5*shapeScaleFactor;
options.bodyAxes.magY = 2.0*shapeScaleFactor;
options.bodyAxes.magZ = 2.0*shapeScaleFactor;
options.bodyAxes.lineWidth = 2.5;

% helper lines
options.helperLines.show = true;
options.helperLines.lineStyle = ':';
options.helperLines.lineColor = 'k';
options.helperLines.lineWidth = 1.5;

% trajectory
options.trajectory.show = true;
options.trajectory.lineStyle = '-';
options.trajectory.lineColor = 'k';
options.trajectory.lineWidth = 1.5;

%% Plot body and trajectory
plotTrajectoryAndBodyE(h_fig3, shape, mXYZe, mEulerAngles, options);

%% Plot Earth axes
% hold on;
% xMax = max([max(abs(mXYZe(:,1))),5]);
% yMax = max([max(abs(mXYZe(:,2))),5]);
% zMax = 0.05*xMax; % max([abs(max(vXYZe(1))),0.18*xMax]);
% vXYZ0 = [0,0,0];
% vExtent = [xMax,yMax,zMax];
% plotEarthAxes(h_fig3, vXYZ0, vExtent);
% xlabel('x_E (m)'); ylabel('y_E (m)'); zlabel('z_E (m)')
% hold off