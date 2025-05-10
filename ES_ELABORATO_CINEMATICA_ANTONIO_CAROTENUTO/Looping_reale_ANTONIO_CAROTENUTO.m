clear all; close all; clc;
%% condizioni iniziali
phi0=0; 
psi0=0; 
theta0=0;
X0=0;
Y0=0;
Z0=0;
p_max=0;%rad/s
q_max=1;%rad/s
r_max=0;%rad/s
t_fin=5/2*pi;
N=100;
t_span=linspace(0,t_fin,N);
%% Leggi temporali delle componenti della velocita ' angolare del velivolo nel riferimento body
p=@(t) 0;

vBreakPointsQ(1,:)=[0, t_fin/5, t_fin*4/5, t_fin]; %vettore dei tempi di stop
vBreakPointsQ(2,:)=[0, q_max, q_max,0]; %vettore del valore della velocità q 

q=@(t)...
interp1(vBreakPointsQ(1 ,:) ,vBreakPointsQ(2 ,:) ,t);
q_plot=q(linspace(0,t_fin,N));
r=@(t) 0;

%% Leggi temporali delle componenti delle velocità di traslazione nel riferimento body
u0=100;%m/s

vBreakPointsU(1,:)=[0, t_fin/5, t_fin*4/5, t_fin ]; %vettore dei tempi di stop
vBreakPointsU(2 ,:)=[u0+10, u0, u0,u0];%vettore del valore della velocità u 


u=@(t)...
interp1(vBreakPointsU(1 ,:) ,vBreakPointsU(2 ,:) ,t);

v0=0;%m/s
w0 = 0 ; %m/s

v=@(t) 0;
w=@(t) 0;
%% risoluzione della Gimbal Equation
dQuatdt=@(t ,Q)...
0.5 * [ 0 -p(t) -q(t) -r(t) ;
p(t) 0 r(t) q(t) ; 
q(t) -r(t) 0 p(t);
r(t) q(t) -p(t) 0]*Q;

options=odeset ('RelTol' ,1e-9,'AbsTol' ,1e-9*ones (1,4) ) ; 
Q0=angle2quat ( psi0 , theta0 , phi0 ) ;
[vTime,vQuat]=ode45(dQuatdt,[0 t_fin],Q0,options);

%% Grafica per l 'orientamento
quat0 = vQuat(:,1);
quatx = vQuat(:,2);
quaty = vQuat(:,3);
quatz = vQuat(:,4);
figure (1) ;
subplot 121
plot (vTime , quat0 , 'k.-' , vTime , quatx , 'b.-' , vTime , quaty , 'c-',vTime, quatz,'r.-');
legend('q 0', 'q x', 'q y','q z'); title('Quaternion components');
xlabel ('t') ; ylabel ('Quaternions') ;

[ vpsi , vtheta , vphi ] =quat2angle (vQuat) ;
subplot 122
plot(vTime, vpsi,'k.-', vTime, vtheta, 'b.-', vTime, vphi,'r-');
legend('\psi', '\theta', '\phi');
title('Euler angles'); xlabel('t'); ylabel('Euler angles');
%% risoluzione della Navigation Equation
%Funzione interpolante per conoscere la storia del quaternione
Quat=@(t)...
[interp1( vTime , vQuat (:,1),t),...
interp1(vTime,vQuat(:,2),t) ,... 
interp1(vTime,vQuat(:,3),t) ,... 
interp1(vTime,vQuat(:,4),t) ];

%Matrice per passare dagli assi terra agli assi body
T_BE=@(Q)... 
    quat2dcm (Q) ;%quat2dcm calcola la matrice dei coseni direttori data una matrice
%di quaternioni

dPosEdt=@(t,PosE)... 
transpose(quat2dcm(Quat(t)))*[u(t);v(t);w(t) ];%il primo pezzo è TEB
options=odeset ('RelTol' ,1e-9,'AbsTol' ,1e-9*ones (3 ,1) ) ;
PosE0=[0; 0; 0];
[vTime2,vPosE]=ode45(dPosEdt,vTime,PosE0,options);
vXe=vPosE(: ,1) ; 
vYe=vPosE(: ,2) ;
vZe=vPosE(:,3)+Z0;
%% Figura
figure
plot(t_span,q_plot,'b');
%% Setup the figure/scene for 3D visualization
h_fig3 = figure (3) ;
grid on;
hold on;
light('Position' ,[1 0 -4],'Style','local');
% Trick to have Ze pointing downward and correct visualization
set(gca,'XDir','reverse'); set(gca,'ZDir','reverse'); daspect ([1 1 1]) ;
%% Load aircraft shape
shapeScaleFactor = 50.0;
%shape = loadAircraftMAT('aircraft pa24 =250.mat', scale factor); 
shape = loadAircraftMAT('aircraft_mig29.mat', shapeScaleFactor);
mXYZe = [vPosE(: ,1) ,vPosE(: ,2) ,vPosE(: ,3)+Z0]; mEulerAngles = [ vpsi , vtheta , vphi ] ;
%% Settings
% General settings
options.samples = [1, 61,111,141:50:numel(vTime)]; % [1,40,80,120,160,200,250,300,320,numel(vTime) ];
options.theView = [105 15];

% body axes settings
options.bodyAxes.show = true;
options . bodyAxes .magX = 1* shapeScaleFactor ;
options . bodyAxes .magY = 1* shapeScaleFactor ; options . bodyAxes .magZ = 1* shapeScaleFactor ; options.bodyAxes.lineWidth = 1.5;
% helper lines
options . helperLines .show = true ;
options . helperLines . lineStyle = ':';
options . helperLines . lineColor = 'w' ;
options . helperLines . lineWidth = 0.00001;
% trajectory
options . trajectory .show = true ; 
options . trajectory . lineStyle = '-';
options . trajectory . lineColor = 'k'; options . trajectory . lineWidth = 1.5;
%% Plot body and trajector
plotTrajectoryAndBodyE (h_fig3,shape,mXYZe,mEulerAngles,options);
%% Plot Earth axes
hold on; xMax=max([max(abs(mXYZe(:,1))),5]); yMax=max([max(abs(mXYZe(:,2))),5]); zMax = 0.05*xMax;
vXYZ0= [0,0,0];
vExtent = [xMax,yMax,zMax];
plotEarthAxes(h_fig3 , vXYZ0, vExtent);
xlabel ('x E (m)') ; ylabel ('y E (m)') ;
zlabel ('z E (m)') 
hold off
axis([-50 280 -100 100 -250 40]);