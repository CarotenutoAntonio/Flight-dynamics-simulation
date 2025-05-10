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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODIFICO QUA %%

%% VOGLIO FARE STUDIO PARAMETRICO, VEDO COSA ACCADE AL VARIARE DELLA POS
%% DEL BARICENTRO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %  Allocazione in memoria delle matrici di raccolta dati
    trim_matrix = zeros(5,4);
    
    X_CG_vector = [myAC.Xcg_adim-0.05,myAC.Xcg_adim,myAC.Xcg_adim+0.05,0.42];
    for i = 1:4
        
        myAC.Xcg_adim = X_CG_vector(i);
        myAC.Cm_alpha = -myAC.CL_alpha*(myAC.Xn_adim - myAC.Xcg_adim);


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

    t_1=1; %istante in cui si lasciano i comandi 
    t_fin=30; %istante finale 
    state0 = [V0,alpha0_rad,q0,xEG_0,zEG_0,theta0_rad]; %vettore di stato iniziale
    %della dinamica a comandi liberi, ha come componenti gli elementi della
    %condizione di trim 

    % Assegnazione delle leggi temporali dei comandi di volo    
    delta_e = @(t) interp1([0,  t_1],...
                           [delta_e0_rad, delta_e0_rad],...
                           t,'linear');
    delta_s = @(t) interp1([0,t_fin],[delta_s0_rad,delta_s0_rad],t,'linear');%cost
    delta_T = @(t) interp1([0,t_fin],[delta_T0,delta_T0],t,'linear');%cost
    delta_tab = @(t) 0*t; %fissato a 0 
%disp('')
%disp('Condizione di trim:')
%disp(['Velocità V_0= ',num2str(V0), 'm/s'])
%disp(['Angolo d%attacco alpha_0=' ,num2str(alpha_0_deg),'deg'])
%disp(['Elevatore delta_e_0= ',num2str(delta_e_0_deg),'deg'])
%disp(['Stabilizzatore delta_s_0= ',num2str(delta_s_0_deg),'deg'])
%disp(['Manetta delta_T_0= ',num2str(delta_T0)])
%disp(['Angolo del Tab=0',0])

%% PASSO2 Integrazione dell’istante t0=0s all’istante t1=1s comandi bloccati
%Integrazione a comandi bloccati 
[T1,Y1]=ode45(@eqLongDynamicStickFixed_,[0 t_1],state0);
%estrazione delle variabili
V1=Y1(:,1);
alpha1=Y1(:,2);
q1=Y1(:,3);
xeg1=Y1(:,4);
zeg1=Y1(:,5);
gamma1=Y1(:,6);
theta1=gamma1+alpha1-myAC.mu_x;
% Tutti i vettori ottenuti sono le condizioni di trim protratte fino all'istante 1 
% prevedibilemente sono costanti 
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


%% modifico qua per costruire matrice m_state
m_state=[T,V,alpha,q,xeg,zeg,theta,gamma,dot_delta_e,delta_e_];

%% creato matrice di stato da t=0 a t=tfin


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


%% modifico qua
    %% costruisco una matriceper le altre variabili di interesse
    %% fattori di carico e coefficiente di momento e momento
    m_exit=[T,fza,fxa,v_Ch_e,v_H_Ae,q_dot_deg_s];

%%
    if i == 1
           A_state = m_state;
           A_exit=m_exit;

        elseif i == 2
               B_state = m_state;
               B_exit=m_exit;

        elseif i == 3
               C_state = m_state;
               C_exit= m_exit;

        elseif i == 4
               D_state = m_state;
               D_exit=m_exit;
               
     end

    end
    


 end



%% GRAFICA
% VARIABILI DI STATO
figure(1)
subplot 311 % X_EG
plot(A_state(:,1),A_state(:,5),'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5);
hold on;
plot(B_state(:,1),B_state(:,5),'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
plot(C_state(:,1),C_state(:,5),':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
plot(D_state(:,1),D_state(:,5),'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',1.5);
grid on
legend('Location','northwest')
lgd = legend('$\hat{X}_{CG} = 0.27$','$\hat{X}_{CG} = 0.32$','$\hat{X}_{CG} = 0.37$','$\hat{X}_{CG} = 0.43$');
lgd.Interpreter = 'latex'; 
xlim([0 t_fin])
xlabel('t (s)','fontsize',12);
ylim([0 6000])
ylabel('x_{EG} (m)','fontsize',12);


subplot 312 % VARIAZIONE DI QUOTA
plot(A_state(:,1),A_state(:,6)-zEG_0,'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5);
hold on;
plot(B_state(:,1),B_state(:,6)-zEG_0,'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
plot(C_state(:,1),C_state(:,6)-zEG_0,':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
plot(D_state(:,1),D_state(:,6)-zEG_0,'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',1.5);
grid on
legend('Location','northwest')
lgd = legend('$\hat{X}_{CG} = 0.27$','$\hat{X}_{CG} = 0.32$','$\hat{X}_{CG} = 0.37$','$\hat{X}_{CG} = 0.43$');
lgd.Interpreter = 'latex'; 
xlim([0 t_fin])
xlabel('t (s)','fontsize',12);
ylim([0 6000])
ylabel('\Delta h (m)','fontsize',12)


subplot 313 % VARIAZIONE DI VELOCITA'
plot(A_state(:,1),A_state(:,2)-V0,'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5);
hold on;
plot(B_state(:,1),B_state(:,2)-V0,'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
plot(C_state(:,1),C_state(:,2)-V0,':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
plot(D_state(:,1),D_state(:,2)-V0,'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',1.5);
grid on
legend('Location','southeast')
lgd = legend('$\hat{X}_{CG} = 0.27$','$\hat{X}_{CG} = 0.32$','$\hat{X}_{CG} = 0.37$','$\hat{X}_{CG} = 0.42$');
lgd.Interpreter = 'latex'; 
xlim([0 t_fin])
xlabel('t (s)','fontsize',12);
ylim([-200 100])
ylabel('\Delta V (m/s)','fontsize',12)



% VARIABILI DI STATO
figure(2)
subplot 221 % VARIAZIONE DI ALPHA
plot(A_state(:,1),convang(A_state(:,3)-A_state(1,3),'rad','deg'),'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5);
hold on;
plot(B_state(:,1),convang(B_state(:,3)-B_state(1,3),'rad','deg'),'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
plot(C_state(:,1),convang(C_state(:,3)-C_state(1,3),'rad','deg'),':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
plot(D_state(:,1),convang(D_state(:,3)-D_state(1,3),'rad','deg'),'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',1.5);
grid on
%legend('Location','northwest')
lgd = legend('$\hat{X}_{CG} = 0.27$','$\hat{X}_{CG} = 0.32$','$\hat{X}_{CG} = 0.37$','$\hat{X}_{CG} = 0.43$');
lgd.Interpreter = 'latex'; 
xlim([0 6])
ylim([-10 5]) 
ylabel('\Delta \alpha (deg)','fontsize',12);


subplot 222 % VARIAZIONE DI ALPHA
plot(A_state(:,1),convang(A_state(:,3)-A_state(1,3),'rad','deg'),'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5);
hold on;
plot(B_state(:,1),convang(B_state(:,3)-B_state(1,3),'rad','deg'),'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
plot(C_state(:,1),convang(C_state(:,3)-C_state(1,3),'rad','deg'),':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
plot(D_state(:,1),convang(D_state(:,3)-D_state(1,3),'rad','deg'),'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',1.5);
grid on
%legend('Location','northwest')
lgd = legend('$\hat{X}_{CG} = 0.27$','$\hat{X}_{CG} = 0.32$','$\hat{X}_{CG} = 0.37$','$\hat{X}_{CG} = 0.43$');
lgd.Interpreter = 'latex'; 
xlim([0 t_fin])
ylim([-10 5]) 
ylabel('\Delta \alpha (deg)','fontsize',12);


subplot 223 % VARIAZIONE DI Q (VELOCITà ANGOLARE)
plot(A_state(:,1),convangvel(A_state(:,4)-A_state(1,4),'rad/s','deg/s'),'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5);
hold on;
plot(B_state(:,1),convangvel(B_state(:,4)-B_state(1,4),'rad/s','deg/s'),'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
plot(C_state(:,1),convangvel(C_state(:,4)-C_state(1,4),'rad/s','deg/s'),':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
plot(D_state(:,1),convangvel(D_state(:,4)-D_state(1,4),'rad/s','deg/s'),'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',1.5);
grid on
%legend('Location','northwest')
lgd = legend('$\hat{X}_{CG} = 0.27$','$\hat{X}_{CG} = 0.32$','$\hat{X}_{CG} = 0.37$','$\hat{X}_{CG} = 0.43$');
lgd.Interpreter = 'latex'; 
xlim([0 5])
xlabel('t (s)','fontsize',12);
ylim([-20 20]) 
ylabel('\Delta q (deg/s)','fontsize',12);

subplot 224 % VARIAZIONE DI VELOCITà ANGOLARE
plot(A_state(:,1),convangvel(A_state(:,4)-A_state(1,4),'rad/s','deg/s'),'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5);
hold on;
plot(B_state(:,1),convangvel(B_state(:,4)-B_state(1,4),'rad/s','deg/s'),'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
plot(C_state(:,1),convangvel(C_state(:,4)-C_state(1,4),'rad/s','deg/s'),':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
plot(D_state(:,1),convangvel(D_state(:,4)-D_state(1,4),'rad/s','deg/s'),'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',1.5);
grid on
%legend('Location','northwest')
lgd = legend('$\hat{X}_{CG} = 0.27$','$\hat{X}_{CG} = 0.32$','$\hat{X}_{CG} = 0.37$','$\hat{X}_{CG} = 0.43$');
lgd.Interpreter = 'latex';
xlim([0 t_fin])
xlabel('t (s)','fontsize',12);
ylim([-20 10]) 
ylabel('\Delta q (deg/s)','fontsize',12);


figure(3)
% modifico qua abbastanza%
subplot 221
plot(A_state(:,1),convang(A_state(:,7),'rad','deg'),'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
hold on;
plot(A_state(:,1),convang(A_state(:,8),'rad','deg'),':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
grid on
lgd = legend('$\theta(t)$','$\gamma(t)$');
lgd.Interpreter = 'latex'; 
lgd.FontSize = 11;
xlim([0 5])
xlabel('t (s)','fontsize',12);
ylim([-10 10])
ylabel('$(deg)$','interpreter','latex','fontsize',11)
title('$\hat{X}_{CG} = 0.27$','interpreter','latex','fontsize',11)

subplot 222
plot(B_state(:,1),convang(B_state(:,7),'rad','deg'),'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
hold on;
plot(B_state(:,1),convang(B_state(:,8),'rad','deg'),':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
grid on
lgd = legend('$\theta(t)$','$\gamma(t)$');
lgd.Interpreter = 'latex'; 
lgd.FontSize = 11;
xlim([0 5])
xlabel('t (s)','fontsize',12);
ylim([-10 10])
ylabel('$(deg)$','interpreter','latex','fontsize',11)
title('$\hat{X}_{CG} = 0.32$','interpreter','latex','fontsize',11)


subplot 223
plot(C_state(:,1),convang(C_state(:,7),'rad','deg'),'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
hold on;
plot(C_state(:,1),convang(C_state(:,8),'rad','deg'),':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
grid on
lgd = legend('$\theta(t)$','$\gamma(t)$');
lgd.Interpreter = 'latex'; 
lgd.FontSize = 11;
xlim([0 5])
xlabel('t (s)','fontsize',12);
ylim([-10 10])
ylabel('$(deg)$','interpreter','latex','fontsize',11)
title('$\hat{X}_{CG} = 0.37$','interpreter','latex','fontsize',11)

subplot 224
plot(D_state(:,1),convang(D_state(:,7),'rad','deg'),'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
hold on;
plot(D_state(:,1),convang(D_state(:,8),'rad','deg'),':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
grid on
lgd = legend('$\theta(t)$','$\gamma(t)$');
lgd.Interpreter = 'latex'; 
lgd.FontSize = 11;
xlim([0 5])
xlabel('t (s)','fontsize',12);
ylim([-10 10])
ylabel('$(deg)$','interpreter','latex','fontsize',11)
title('$\hat{X}_{CG} = 0.43$','interpreter','latex','fontsize',11)




% PARAMETRI PER LO STUDIO DELLA CERNIERA
figure(4)
subplot 211 % DELTA E
plot(A_state(:,1),convang(A_state(:,10),'rad','deg'),'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5);
hold on;
plot(B_state(:,1),convang(B_state(:,10),'rad','deg'),'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
plot(C_state(:,1),convang(C_state(:,9),'rad','deg'),':','color',[0.4660, 0.6740, 0.1880],'LineWidth',0.5);
plot(D_state(:,1),convang(D_state(:,9),'rad','deg'),'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',0.5);
grid on
legend('Location','southeast')
lgd = legend('$\hat{X}_{CG} = 0.27$','$\hat{X}_{CG} = 0.32$','$\hat{X}_{CG} = 0.37$','$\hat{X}_{CG} = 0.43$');
lgd.Interpreter = 'latex'; 
xlim([0 5])
ylim([-10 10]) 
xlabel('t (s)','fontsize',12);
ylabel('\delta_{e} (deg)','fontsize',12);


subplot 212 % DELTA E DOT
plot(A_state(:,1),convang(A_state(:,9),'rad','deg'),'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5);
hold on;
plot(B_state(:,1),convang(B_state(:,9),'rad','deg'),'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
plot(C_state(:,1),convang(C_state(:,9),'rad','deg'),':','color',[0.4660, 0.6740, 0.1880],'LineWidth',0.5);
plot(D_state(:,1),convang(D_state(:,9),'rad','deg'),'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',0.5);
grid on
legend('Location','southeast')
lgd = legend('$\hat{X}_{CG} = 0.27$','$\hat{X}_{CG} = 0.32$','$\hat{X}_{CG} = 0.37$','$\hat{X}_{CG} = 0.43$');
lgd.Interpreter = 'latex'; 
xlim([0 5])
ylim([-400 400]) 
ylabel('dot{delta}_{e} (deg/s)','fontsize',12);
xlabel('t (s)','fontsize',12);

figure(5)
subplot 211 % COEFFICIENTE DEL MOMENTO DI CERNIERA
plot(A_exit(:,1),A_exit(:,4),'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5);
hold on;
plot(B_exit(:,1),B_exit(:,4),'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
plot(C_exit(:,1),C_exit(:,4),':','color',[0.4660, 0.6740, 0.1880],'LineWidth',0.5);
plot(D_exit(:,1),D_exit(:,4),'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',0.5);
grid on
legend('Location','southeast')
lgd = legend('$\hat{X}_{CG} = 0.27$','$\hat{X}_{CG} = 0.32$','$\hat{X}_{CG} = 0.37$','$\hat{X}_{CG} = 0.43$');
lgd.Interpreter = 'latex'; 
xlim([0 3])
ylim([-0.04 0.04]) 
ylabel('C_{He}','fontsize',12);
xlabel('t (s)','fontsize',12);

subplot 212 % MOMENTO DI CERNIERA
plot(A_exit(:,1),A_exit(:,5),'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5);
hold on;
plot(B_exit(:,1),B_exit(:,5),'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
plot(C_exit(:,1),C_exit(:,5),':','color',[0.4660, 0.6740, 0.1880],'LineWidth',0.5);
plot(D_exit(:,1),D_exit(:,5),'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',0.5);
grid on
legend('Location','southeast')
lgd = legend('$\hat{X}_{CG} = 0.27$','$\hat{X}_{CG} = 0.32$','$\hat{X}_{CG} = 0.37$','$\hat{X}_{CG} = 0.43$');
lgd.Interpreter = 'latex'; 
xlim([0 3])
xlabel('t (s)','fontsize',12);
ylim([-400 500]) 
ylabel('H_{e} (N m)','fontsize',12);
xlabel('t (s)','fontsize',12);

% FATTORI DI CARICO E ACCELERAZONE ANGOLARE
figure(6);
subplot 221  % FATTORE DI CARICO NORMALE
plot(A_exit(:,1),A_exit(:,2),'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5);
hold on;
plot(B_exit(:,1),B_exit(:,2),'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
plot(C_exit(:,1),C_exit(:,2),':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
plot(D_exit(:,1),D_exit(:,2),'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',1.5);
grid on
legend('Location','northeast')
lgd = legend('$\hat{X}_{CG} = 0.27$','$\hat{X}_{CG} = 0.32$','$\hat{X}_{CG} = 0.37$','$\hat{X}_{CG} = 0.43$');
lgd.Interpreter = 'latex'; 
ylabel('f_{z_A}','fontsize',12);
xlabel('t (s)','fontsize',12);
xlim([0 5]);
ylim([-5 5]);
grid on ;

subplot 222
plot(A_exit(:,1),A_exit(:,2),'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5);
hold on;
plot(B_exit(:,1),B_exit(:,2),'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
plot(C_exit(:,1),C_exit(:,2),':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
plot(D_exit(:,1),D_exit(:,2),'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',1.5);
grid on
legend('Location','northeast')
lgd = legend('$\hat{X}_{CG} = 0.27$','$\hat{X}_{CG} = 0.32$','$\hat{X}_{CG} = 0.37$','$\hat{X}_{CG} = 0.43$');
lgd.Interpreter = 'latex'; 
ylabel('f_{z_A}','fontsize',12);
xlabel('t (s)','fontsize',12);
xlim([0 t_fin]);
ylim([-5 5]);
grid on ;


subplot 223 % FATTORE DI CARICO LUNGO X
plot(A_exit(:,1),A_exit(:,3),'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5);
hold on;
plot(B_exit(:,1),B_exit(:,3),'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
plot(C_exit(:,1),C_exit(:,3),':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
plot(D_exit(:,1),D_exit(:,3),'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',1.5);
grid on
legend('Location','northwest')
lgd = legend('$\hat{X}_{CG} = 0.27$','$\hat{X}_{CG} = 0.32$','$\hat{X}_{CG} = 0.37$','$\hat{X}_{CG} = 0.43$');
lgd.Interpreter = 'latex'; 
ylabel('f_{x_A}','fontsize',12);
xlabel('t (s)','FontSize',12);
xlim([0 5]);
ylim([-1 5])
grid on;

subplot 224 
plot(A_exit(:,1),A_exit(:,3),'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5);
hold on;
plot(B_exit(:,1),B_exit(:,3),'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
plot(C_exit(:,1),C_exit(:,3),':','color',[0.4660, 0.6740, 0.1880],'LineWidth',0.3);
plot(D_exit(:,1),D_exit(:,3),'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',0.3);
grid on
legend('Location','northeast')
lgd = legend('$\hat{X}_{CG} = 0.27$','$\hat{X}_{CG} = 0.32$','$\hat{X}_{CG} = 0.37$','$\hat{X}_{CG} = 0.43$');
lgd.Interpreter = 'latex'; 
ylabel('f_{x_A}','fontsize',12);
xlabel('t (s)','fontsize',12);
xlim([0 t_fin]);
ylim([-1 5])
grid on;

figure(7)
subplot 211 % Q DOT (ACCELEREZIONE ANGOLARE)
plot(A_exit(:,1),A_exit(:,6),'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5);
hold on;
plot(B_exit(:,1),B_exit(:,6),'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
plot(C_exit(:,1),C_exit(:,6),':','color',[0.4660, 0.6740, 0.1880],'LineWidth',0.3);
plot(D_exit(:,1),D_exit(:,6),'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',0.3);
grid on
legend('Location','southeast')
lgd = legend('$\hat{X}_{CG} = 0.27$','$\hat{X}_{CG} = 0.32$','$\hat{X}_{CG} = 0.37$','$\hat{X}_{CG} = 0.43$');
lgd.Interpreter = 'latex'; 
xlabel('t (s)','fontsize',12);
ylabel('dot{q} [deg/s^2]','fontsize',12);
xlim([0 5]);
ylim([-180 50]);
grid on;

subplot 212
plot(A_exit(:,1),A_exit(:,6),'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',2.5);
hold on;
plot(B_exit(:,1),B_exit(:,6),'--','color',[0, 0.4470, 0.7410],'LineWidth',2.5);
plot(C_exit(:,1),C_exit(:,6),':','color',[0.4660, 0.6740, 0.1880],'LineWidth',0.8);
plot(D_exit(:,1),D_exit(:,6),'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',0.8);
grid on
legend('Location','southeast')
lgd = legend('$\hat{X}_{CG} = 0.27$','$\hat{X}_{CG} = 0.32$','$\hat{X}_{CG} = 0.37$','$\hat{X}_{CG} = 0.43$');
lgd.Interpreter = 'latex'; 
xlabel('t (s)','fontsize',12);
ylabel('dot{q} [deg/s^2]','fontsize',12);
xlim([0 t_fin]);
ylim([-180 50]);
grid on;
% 