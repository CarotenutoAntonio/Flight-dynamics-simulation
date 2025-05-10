clc;clear;close all;
clear ; clc; close all;
%% Script di risoluzione del problema di valori iniziali , moto a 3 DoF, comandi bloccati
%% Declarations
global g...                   %Accelerazione di gravita'
    zEG_0 V0 q0 gamma0...  %Condizioni iniziali
    rho0 ...               %Densità dell'aria all'altitudine h = (-zEG_0)
    myAC   ...             %Oggetto 'Velivolo'
    delta_e...
    delta_s...
    delta_T ...
    delta_tab
%% Populate aircraft data
aircraftDataFileName = 'DSV_Aircraft_data.txt';

%% Aircraft object
myAC = DSVAircraft(aircraftDataFileName);

if (myAC.err == -1)
    disp('... terminating.');
else
    disp(['File ', aircraftDataFileName, ' letto correttamente.']);
end
%% Condizioni iniziali
xEG_0=0;
zEG_0 = -4000; % =quota (m)
V0 = 257; % velocita' di volo
q0 = 0; % velocita ' angolare di beccheggio (rad/s)
gamma0 = convang(0, 'deg', 'rad'); % angolo di salita (rad)

%% Densita ' con modello ISA
[air_Temp0, sound_speed0, air_pressure0, rho0] = atmosisa(-zEG_0);

%% Accelerazione di gravita '
g = 9.81; % (m/sˆ2)

%% Tentativo iniziale
x0 = [
    0;  %alpha0
    0;  %deltae0
    0;  %deltas0
    0.5 % delta_T_0
    ];
%% Minimo della funzione di costo
% Aeq, in Aeq*x=beq linear constraint
Aeq = zeros (4, 4);
ceq = zeros (4, 1);
Aeq(3, 3) = 1;
delta_s_0 = convang(-0.0, 'deg', 'rad');
ceq(3, 1) = delta_s_0; %tiene deleta_s fissato
% bounds
lb = [convang(-15, 'deg', 'rad'), ... % minimo alpha
    convang(-20, 'deg', 'rad'), ... % minima deflessione dell 'equilibratore
    convang(-5, 'deg', 'rad'), ... % minima incidenza dello stabilizzatore
    0.2]; ... % minima manetta
    ub = [convang(15, 'deg', 'rad'), ... % massimo alpha
    convang(13, 'deg', 'rad'), ... % massima deflessione dell 'equilibratore
    convang(2, 'deg', 'rad'), ... % massima incidenza dello stabilizzatore
    1.0]; %massima manetta
options = optimset( ...
    'tolfun', 1e-9, ... %treshold
    'Algorithm', 'interior-point' ... % algor. type
    );
[x, fval] = ...
    fmincon(@costLongEquilibriumStaticStickFixed, ...
    x0, ...
    [], ... %A,A*x<=b
    [], ... %b
    Aeq, ... % Aeq , Aeq*x=beq
    ceq, ... % beq
    lb, ub, ...
    @myNonLinearConstraint, ...
    options);
alpha0_rad = x(1);
alpha_0_deg = convang(alpha0_rad, 'rad', 'deg');

theta0_rad = alpha0_rad - myAC.mu_x + gamma0;
theta_0_deg = convang(theta0_rad, 'rad', 'deg');

delta_e0_rad = x(2);
delta_e_0_deg = convang(delta_e0_rad, 'rad', 'deg');

delta_s0_rad = x(3);
delta_s_0_deg = convang(x(3), 'rad', 'deg');
delta_T0 = x(4);

% Calcolo Trim Tab 
theta_trim = gamma0 + alpha0_rad - myAC.mu_x;
alpha_body_0 = alpha0_rad-myAC.mu_x;               % alpha body
alpha_H_0 = (1-myAC.DepsDalpha)*(alpha_body_0)...   % alpha H
             -myAC.eps_0 + delta_s0_rad + myAC.mu_x;
CH_A_e_delta_tab_segnato = ((myAC.mass_e*myAC.ec_adim*myAC.mac_e)...
         *g*cos(theta_trim))/(0.5*rho0*V0^2*myAC.S_e*myAC.mac_e)...
         - (myAC.Ch_e_0 + myAC.Ch_e_alpha*alpha_H_0 ...
         + myAC.Ch_e_delta_e*delta_e0_rad ...
         + myAC.Ch_e_delta_s*delta_s0_rad);
delta_tab_segnato = CH_A_e_delta_tab_segnato/myAC.Ch_e_delta_tab;
%% modifico qua faccio due leggi della manetta
for i=1:3
    t_1=1;
    t_2=10;
    t_fin=150;

    if i==1

    t_3=15;
    t_4=20;
    delta_tab = @(t) interp1([0, t_1, t_2, t_3, t_4, t_fin], ...
        [0, delta_tab_segnato, delta_tab_segnato, delta_tab_segnato*3/5, delta_tab_segnato, delta_tab_segnato] ...
        , t,'linear');
    elseif i==2

    t_3=20; %25
    t_4=30; %40
    delta_tab = @(t) interp1([0, t_1, t_2, t_3, t_4, t_fin], ...
        [0, delta_tab_segnato, delta_tab_segnato, delta_tab_segnato*3/5, delta_tab_segnato, delta_tab_segnato] ...
        , t,'linear');
    elseif i==3
    t_3=15; %t_3=20
    t_4=20;
    delta_tab = @(t) interp1([0, t_1, t_2, t_3, t_4, t_fin], ...
        [0, delta_tab_segnato, delta_tab_segnato, delta_tab_segnato*1/5, delta_tab_segnato, delta_tab_segnato] ...
        , t,'linear');
    end
   

t_1=1;
t_2=10;
t_fin=150;
state0 = [V0,alpha0_rad,q0,xEG_0,zEG_0,theta0_rad];
% Assegnazione delle leggi temporali dei comandi di volo
delta_e = @(t) interp1([0,  t_1],...
    [delta_e0_rad, delta_e0_rad],...
    t,'linear');
delta_s = @(t) interp1([0,t_fin],[delta_s0_rad,delta_s0_rad],t,'linear');
delta_T = @(t) interp1([0,t_fin],[delta_T0,delta_T0],t,'linear');
disp('')
disp('Condizione di trim:')
disp(['Velocità V_0= ',num2str(V0), 'm/s'])
disp(['Angolo d’attacco alpha_0=' ,num2str(alpha_0_deg),'deg'])
disp(['Elevatore delta_e_0= ',num2str(delta_e_0_deg),'deg'])
disp(['Stabilizzatore delta_s_0= ',num2str(delta_s_0_deg),'deg'])
disp(['Manetta delta_T_0= ',num2str(delta_T0)])
disp(['Angolo del Tab=0 ',0])

%% PASSO2 Integrazione dell’istante t0=0s all’istante t1=1s comandi bloccati

[T1,Y1]=ode45(@eqLongDynamicStickFixed_,[0 t_1],state0);
V1=Y1(:,1);
alpha1=Y1(:,2);
q1=Y1(:,3);
xeg1=Y1(:,4);
zeg1=Y1(:,5);
gamma1=Y1(:,6);
theta1=gamma1+alpha1-myAC.mu_x;
%% PASSO3 Risoluzione del sistema di equazioni 7.42, dall’istante t1=1s a
% tf=30s per il caso di comandi liberi

delta_e0_dot_rad = convangvel(0,'deg/s','rad/s');
state_2 = [V1(end),alpha1(end),q1(end),xeg1(end),zeg1(end),theta1(end),...
    delta_e0_dot_rad,delta_e(T1(end))];
[T2,Y2] = ode45(@eqLongDynamicStickFree_,[t_1 t_fin],state_2);

%Estrazione grandezze di interesse
V2=Y2(:,1);
alpha2=Y2(:,2);
q2=Y2(:,3);
xeg2=Y2(:,4);
zeg2=Y2(:,5);
theta2=Y2(:,6);
delta_e2_dot_rad=Y2(:,7);
delta_e2=Y2(:,8);

% Mi aspetto che tutti questi valori siano costanti poichè tramite
% l'introduzione della tab nel primo secondo in cui i comandi osno bloccati
% il velivolo risulta trimmato a comandi liberi 

%Composizione dei vettori di stato di tutta la manovra:(vettori colonna)
T=[T1; T2]; %[s]
V=[V1; V2]; %m/s
alpha=[alpha1; alpha2]; %rad
q=[q1; q2]; %rad/s
xeg=[xeg1; xeg2]; %m
zeg=[zeg1; zeg2]; %m
theta=[theta1; theta2]; %rad
gamma=theta-alpha+myAC.mu_x; %[rad]
dot_delta_e=[zeros(length(T1),1); delta_e2_dot_rad];  %rad/s
delta_e_=[delta_e0_rad*ones(length(T1),1); delta_e2]; %rad
alpha_dot=gradient(alpha,T); %rad/s
v_dot=gradient(V,T);
q_dot=gradient(q,T); %[rad/s]
q_dot_deg_s=convangvel(q_dot,'rad/s','deg/s'); %deg/s
gamma_dot=gradient(gamma,T); %rad/s
fxa = -sin(gamma)-v_dot./g;
fza    = cos(gamma)+(gamma_dot.*V)./g;
% Angolo d'attacco del piano di coda orizzontale
alpha_H_rad = (1 - myAC.DepsDalpha)*(alpha-myAC.mu_x) - myAC.eps_0 +...
    delta_s(T) + myAC.mu_x; %rad
alpha_H_rad_dot = (1 - myAC.DepsDalpha)*(alpha_dot); %rad/s

% Coefficiente momento di Ceniera
v_Ch_e = myAC.Ch_e_0 + myAC.Ch_e_alpha*alpha_H_rad +...
    myAC.Ch_e_delta_e*delta_e_+...
    myAC.Ch_e_delta_s*delta_s(T)+...
    myAC.Ch_e_delta_tab*delta_tab(T)+...
    (myAC.mac_e/(2*V))...
    *(myAC.Ch_e_alpha_dot*alpha_H_rad_dot+...
    myAC.Ch_e_q*q+...
    myAC.Ch_e_delta_e_dot*dot_delta_e);
% Momento aerodinamico di cerniera
v_H_Ae = 0.5*density(-zeg).*V.^2*myAC.S_e*myAC.mac_e.*v_Ch_e;
% %% GRAFICA


%% modifico qua
delta_tab_rad = delta_tab(T);
delta_tab_deg=convang (delta_tab_rad, 'rad', 'deg');

%Deflessione Tab
figure(1)
plot( T, delta_tab_deg);
grid on;
hold on;
ylabel('\delta_{t.0} (deg)');
% VARIABILI DI STATO

figure(2)
subplot 311 % X_EG
plot(T , xeg,'LineWidth',1.5);
grid on;
hold on;
xlim([0 t_fin])
xlabel('t (s)','fontsize',12);
ylim([0 40000])
ylabel('x_{EG} (m)','fontsize',12);


subplot 312 % VARIAZIONE DI QUOTA
plot(T,(zeg - zEG_0),'LineWidth',1.5);
grid on;
hold on;
xlim([0 t_fin])
xlabel('t (s)','fontsize',12);
ylim([-1000 4000])
ylabel('\Delta h (m)','fontsize',12)


subplot 313 % VARIAZIONE DI VELOCITA'
plot(T,(V- V0),'LineWidth',1.5);
grid on;
hold on;
xlim([0 t_fin])
xlabel('t (s)','fontsize',12);
ylim([-75 75])
ylabel('\Delta V (m/s)','fontsize',12)



% VARIABILI DI STATO
figure(3)
subplot 321 % VARIAZIONE DI ALPHA
plot(T,convang((alpha- alpha0_rad),'rad','deg'),'linewidth',1.5);
grid on
hold on
xlim([0 6])
ylim([-6 10])
ylabel('\Delta \alpha (deg)','fontsize',12);


subplot 322 % VARIAZIONE DI ALPHA
plot(T,convang((alpha- alpha0_rad),'rad','deg'),'linewidth',1.5);
grid on
hold on
xlim([0 t_fin])
ylim([-3 3])
ylabel('\Delta \alpha (deg)','fontsize',12);


subplot 323 % VARIAZIONE DI Q (VELOCITà ANGOLARE)
plot(T,convangvel((q- q0),'rad/s','deg/s'),'LineWidth',1.5);
grid on
hold on
xlim([0 5])
xlabel('t (s)','fontsize',12);
ylim([-5 5])
ylabel('\Delta q (deg/s)','fontsize',12);

subplot 324 % VARIAZIONE DI VELOCITà ANGOLARE
plot(T,convangvel((q- q0),'rad/s','deg/s'),'LineWidth',1.5);
grid on
hold on
xlim([0 t_fin])
xlabel('t (s)','fontsize',12);
ylim([-5 5])
ylabel('\Delta q (deg/s)','fontsize',12);


subplot 325 % VARIAZIONE DI THETA 
plot(T,convang((theta- theta0_rad),'rad','deg'),'LineWidth',1.5);

%modifico qua
grid on
hold on
legend('\Delta\theta(t)');
xlim([0 t_fin])
xlabel('t (s)','fontsize',12);
ylim([-35 35])
ylabel('\Delta\theta (deg)','fontsize',12);



subplot 326 % VARIAZIONE DI  GAMMA
plot(T, convang((gamma - gamma0),'rad','deg'),'linewidth',1.5);
grid on
hold on
legend('\Delta\gamma(t)');
xlim([0 t_fin])
ylim([-35 35])
xlabel('t (s)','fontsize',12);
ylabel('\Delta\gamma(deg)','fontsize',12);

% PARAMETRI PER LO STUDIO DELLA CERNIERA
figure(4)
subplot 411 % DELTA E
plot(T,convang(delta_e_,'rad','deg'),'LineWidth',1.5);
grid on
hold on
xlim([0 150])
%ylim([-8 8])
xlabel('t (s)','fontsize',12);
ylabel('\delta_{e} (deg)','fontsize',12);


subplot 412 % DELTA E DOT
plot(T,convangvel(dot_delta_e,'rad/s','deg/s'),'LineWidth',1.5);
grid on
hold on
xlim([0 150])
ylim([-2 2])
ylabel('dot{delta}_{e} (deg/s)','fontsize',12);
xlabel('t (s)','fontsize',12);

subplot 413 % COEFFICIENTE DEL MOMENTO DI CERNIERA
plot(T,v_Ch_e,'LineWidth',1.5);
grid on
hold on
xlim([0 150])
ylim([-0.02 0.04])
ylabel('C_{He}','fontsize',12);
xlabel('t (s)','fontsize',12);

subplot 414 % MOMENTO DI CERNIERA
plot(T,v_H_Ae,'LineWidth',1.5);
grid on
hold on
xlim([0 150])
xlabel('t (s)','fontsize',12);
ylim([-400 500])
ylabel('H_{e} (N m)','fontsize',12);
xlabel('t (s)','fontsize',12);

% FATTORI DI CARICO E ACCELERAZONE ANGOLARE
figure(5);
subplot 321  % FATTORE DI CARICO NORMALE
plot(T,fza,'linewidth',1.2);
ylabel('f_{z_A}','fontsize',12);
xlabel('t (s)','fontsize',12);
xlim([0 5]);
ylim ([-2 2]);
grid on ;
hold on;

subplot 322
plot(T,fza,'linewidth',1.2);
ylabel('f_{z_A}','fontsize',12);
xlabel('t (s)','fontsize',12);
xlim([0 t_fin]);
ylim ([-1 2]);
grid on ;
hold on;


subplot 323 % FATTORE DI CARICO LUNGO X
plot(T,fxa,'linewidth',1.2);
ylabel('f_{x_A}','fontsize',12);
xlabel('t (s)','FontSize',12);
xlim([0 5]);
ylim ([-2 2]);
grid on;
hold on;

subplot 324
plot(T,fxa,'linewidth',1.2);
ylabel('f_{x_A}','fontsize',12);
xlabel('t (s)','fontsize',12);
xlim([0 t_fin]);
ylim ([-1 1]);
grid on;
hold on;


subplot 325 % Q DOT (ACCELEREZIONE ANGOLARE)
plot(T,q_dot_deg_s,'linewidth',1.2);
xlabel('t (s)','fontsize',12);
ylabel('dot{q} [deg/s^2]','fontsize',12);
xlim([0 5]);
ylim([-2 2]);
grid on;
hold on;

subplot 326
plot(T,q_dot_deg_s,'linewidth',1.2);
xlabel('t (s)','fontsize',12);
ylabel('dot{q} [deg/s^2]','fontsize',12);
xlim([0 t_fin]);
ylim([-1 1.8]);
grid on;
hold on;

end