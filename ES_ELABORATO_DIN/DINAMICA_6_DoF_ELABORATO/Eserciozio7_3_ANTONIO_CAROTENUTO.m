%% Simulazione del moto a 6-DOf a partire da condizioni di trim
clear all; close all; clc;

disp('Moto del velivolo a 6 gradi di libertà');
%% Dichiarazione delle variabili globali
global g...                  %Accelerazione di gravità 
       myAC                  %Oggetto 'Velivolo'

%  Definizione della classe DSVAircraft e dell'oggetto 'Velivolo'
aircraftDataFileName = 'HARV.txt';
myAC = DSVAircraft(aircraftDataFileName);

%  Costanti e condizioni iniziali   
g = 9.81; %Accelerazione di gravità [m/s^2]
xEG_0 = 0; %[m]
yEG_0 = 0; %[m]
zEG_0 = 0; %Altitudine [m]
[air_Temp0,sound_speed0,air_pressure0,rho0] = atmosisa(-zEG_0);

%%%%%% ANGOLI DI EULERO E QUATERNIONE INIZIALI (MIO INTERVENTO)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi0 = convang(0.000,'deg','rad'); %Angolo di bank
theta0 = convang(4.460,'deg','rad'); %Angolo di elevazione
psi0 = convang(0.000,'deg','rad'); %Angolo di heading
eul0=[psi0,theta0,phi0]
quat_0=angle2quat(psi0,theta0,psi0);
q1_0=quat_0(1);
q2_0=quat_0(2);
q3_0=quat_0(3);
q4_0=quat_0(4);

% velocità angolari iniziali
p0 = convangvel(0.000,'deg/s','rad/s'); %Velocità angolare di rollio
q0 = convangvel(0.000,'deg/s','rad/s'); %Velocità angolare di beccheggio
r0 = convangvel(0.000,'deg/s','rad/s'); %Velocità angolare di imbardata
M0 = 0.246010; %Numero di Mach
V0 = M0*sound_speed0; %Velocità lineare del baricentro
alpha_B0 = convang(4.557754,'deg','rad'); %Angolo di attacco valutato rispetto l'asse x Body
beta_der0 = convang(0.000,'deg','rad'); %Angolo di derapata
u0 = V0*cos(beta_der0)*cos(alpha_B0); %Componente della velocità lineare del baricentro lungo l'asse x Body
v0 = V0*sin(beta_der0); %Componente della velocità lineare del baricentro lungo l'asse y Body
w0 = V0*cos(beta_der0)*sin(alpha_B0); %Componente della velocità lineare del baricentro lungo l'asse z Body

%  Comandi di volo iniziali
delta_a0 = convang(0.000,'deg','rad');
delta_e0 = convang(-10.614717,'deg','rad');
delta_r0 = convang(0.000,'deg','rad'); 
delta_T0 = 0.474918; 
    
%% Integrazione delle equazioni del moto a 6-DoF
t_fin = 230; %Tempo di simulazione [s] (Manovra di avvicinamento alla pista di atterraggio di un aeroporto)
state_0 = [u0,v0,w0,p0,q0,r0,xEG_0,yEG_0,zEG_0,q1_0,q2_0,q3_0,q4_0]; %%%%%%%%%%%%%(MIO INTERVENTO)%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Assegnazione delle leggi temporali dei comandi di volo
global delta_a delta_e delta_r delta_T

%Possibile manovra di avvicinamento alla pista di atterraggio di un aeroporto
delta_a_exc1 = convang(-5,'deg','rad');
delta_a = @(t) interp1([0, 40, 41, 45, 46, t_fin],...
                       [delta_a0, delta_a0, delta_a_exc1, delta_a_exc1, delta_a0, delta_a0],t,'linear');
                   
delta_e_exc1 = convang(-9.5,'deg','rad');
delta_e_exc2 = convang(-10,'deg','rad');
delta_e = @(t) interp1([0, 20, 80, 180, t_fin],...
                       [delta_e0, delta_e0, delta_e_exc1, delta_e_exc2, delta_e_exc2],t,'linear');

delta_r_exc1 = convang(0,'deg','rad');
delta_r = @(t) interp1([0, 200, 203, 205, 208, t_fin],...
                       [delta_r0, delta_r0, delta_r_exc1, delta_r_exc1, delta_r0, delta_r0],t,'linear');                  
                   
delta_T_exc1 = delta_T0 ;
delta_T = @(t) interp1([0, t_fin],[delta_T0 delta_T_exc1],t,'linear');

%  Integrazione del sistema di equazioni differenziali
options = odeset('RelTol',1e-9,'AbsTol',1e-9*ones(1,13));
[vTime,mState] = ode45(@eqLongDynamicStickFixed_6DoF_uvwpqr_EulerQ7n2_PROVOIO,[0 t_fin],state_0,options);

vVel_u = mState(:,1);
vVel_v = mState(:,2);
vVel_w = mState(:,3);
vVel = sqrt(vVel_u.^2 + vVel_v.^2 + vVel_w.^2);
vAlpha_B = convang(atan(vVel_w./vVel_u),'rad','deg');
vBeta = convang(asin(vVel_v./vVel),'rad','deg');
v_p = convangvel(mState(:,4),'rad/s','deg/s');
v_q = convangvel(mState(:,5),'rad/s','deg/s');
v_r = convangvel(mState(:,6),'rad/s','deg/s');
vXe = mState(:,7);
vYe = mState(:,8);
vZe = mState(:,9);
%%%% modifico l'output sono i quaternioni%%%%%%%%%%%%%%%%5
q1 =  mState(:,10);
q2 =  mState(:,11);
q3 =  mState(:,12);
q4 =  mState(:,13);
quat=[q1,q2,q3,q4];
[vPsi,vTheta,vPhi]=quat2angle(quat) ;

vGamma = convang(asin(cos(vAlpha_B).*cos(vBeta).*sin(vTheta)-...
                sin(vBeta).*sin(vPhi).*cos(vTheta)-...
                sin(vAlpha_B).*cos(vBeta).*cos(vPhi).*cos(vTheta)),'rad','deg');
%mQuat = zeros(length(vTime),4);
%for i = 1:length(vTime)
 %   mQuat(i,:) = angle2quat(mState(i,12),mState(i,11),mState(i,10),'ZYX');
%end 
 
%% Grafica
%  Diagrammi dei comandi di volo
figure(1)
subplot 411
plot(vTime,convang(delta_a(vTime),'rad','deg'),'b-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([-25 5]) 
ylabel('$\delta_{a} (deg)$','interpreter','latex','fontsize',11);
subplot 412
plot(vTime,convang(delta_e(vTime),'rad','deg'),'b-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([-11 -9])
ylabel('$\delta_{e} (deg)$','interpreter','latex','fontsize',11);
subplot 413
plot(vTime,convang(delta_r(vTime),'rad','deg'),'b-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([-3 12])
ylabel('$\delta_{r} (deg)$','interpreter','latex','fontsize',11);
subplot 414
plot(vTime,delta_T(vTime),'b-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([0 1])
xlabel('$t (s)$','interpreter','latex','fontsize',11);
ylabel('$\delta_{T}$','interpreter','latex','fontsize',11);

%Storie temporali delle variabili di stato
% figure(2)
% plot(vTime,vVel_u,'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5);
% hold on;
% plot(vTime,vVel_v,'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
% plot(vTime,vVel_w,':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
% grid on
% lgd = legend('$u(t)$','$v(t)$','$w(t)$');
% lgd.Interpreter = 'latex'; 
% lgd.FontSize = 11;
% xlim([0 t_fin])
% xlabel('$t (s)$','interpreter','latex','fontsize',11);
% ylim([-30 150])
% ylabel('$(m/s)$','interpreter','latex','fontsize',11)

figure(3)
subplot 611
plot(vTime,vVel,'-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([80 130])
ylabel('$V (m/s)$','interpreter','latex','fontsize',11)
subplot 612
plot(vTime,vAlpha_B,'-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([-5 10])
ylabel('$\alpha_{B} (deg)$','interpreter','latex','fontsize',11)
subplot 613
plot(vTime,vBeta,'-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([-10 20])
ylabel('$\beta (deg)$','interpreter','latex','fontsize',11)
subplot 614
plot(vTime,v_p,'-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([-5 10])
ylabel('$p (deg/s)$','interpreter','latex','fontsize',11)
subplot 615
plot(vTime,v_q,'-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([-5 5])
ylabel('$q (deg/s)$','interpreter','latex','fontsize',11)
subplot 616
plot(vTime,v_r,'-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
xlabel('$t (s)$','interpreter','latex','fontsize',11);
ylim([-8 10])
ylabel('$r (deg/s)$','interpreter','latex','fontsize',11)
 
figure(4)
subplot 611
plot(vTime,vXe,'-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([-2000 15000])
ylabel('$x_{EG} (m)$','interpreter','latex','fontsize',11)
subplot 612
plot(vTime,vYe,'-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([-2500 15000])
ylabel('$y_{EG} (m)$','interpreter','latex','fontsize',11)
subplot 613
plot(vTime,-(vZe - zEG_0),'-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([-2000 500])
ylabel('$\Delta h (m)$','interpreter','latex','fontsize',11)
subplot 614
plot(vTime,convang(vPhi,'rad','deg'),'-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([-10 50])
ylabel('$\phi (deg)$','interpreter','latex','fontsize',11)
subplot 615
plot(vTime,convang(vTheta,'rad','deg'),'-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([-20 10])
ylabel('$\theta (deg)$','interpreter','latex','fontsize',11)
subplot 616
plot(vTime,convang(vPsi,'rad','deg'),'-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
xlabel('$t (s)$','interpreter','latex','fontsize',11);
ylim([-200 400])
ylabel('$\psi (deg)$','interpreter','latex','fontsize',11)


%rappresento traiettoria
h_fig = figure(5);
theView = [180, 90];
scale_factor = 0.0028;
%step = [1 530 1080 2553 3860 5000 6000 7000 7660 7850 8380 9160 10000];
step = linspace(1,length(vTime),5);
myplotTrajectoryAndBody(h_fig,vXe,vYe,vZe,mState(1:length(vTime),10:13),scale_factor,step,theView);
hold off;

h_fig = figure(6);
theView = [143, 37];
scale_factor = 0.00102;
%step = [1 530 1080 2553 3860 5000 6000 7000 7600 8380 10000];
step = linspace(1,length(vTime),5);
myplotTrajectoryAndBody(h_fig,vXe,vYe,vZe,mState(1:length(vTime),10:13),scale_factor,step,theView);
hold off;