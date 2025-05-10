%% Simulazione del moto a 3-DOf a partire da condizioni di trim
clear all; close all; clc;

disp('Moto del velivolo a 3 gradi di libertà');
disp('Risoluzione del problema di trim ad una data altitudine e velocità di volo');

%% Dichiarazione delle variabili globali
global g...                  %Accelerazione di gravità 
       zEG_0 V0 q0 gamma0... %Condizioni iniziali
       rho0 ...              %Densità dell'aria all'altitudine h = (-zEG_0)
       myAC                  %Oggetto 'Velivolo'

%% Ricerca delle condizioni di trim
aircraftDataFileName = 'DSV_Aircraft_data.txt';

%  Definizione dell'oggetto 'Velivolo'
myAC = DSVAircraft(aircraftDataFileName);

if (myAC.err == -1)
    disp('Terminazione.')
else
    disp(['File ',aircraftDataFileName,' letto correttamente.']);
    
    %  Allocazione in memoria delle matrici di raccolta dati
    trim_matrix = zeros(5,4);
    
    X_CG_vector = [myAC.Xcg_adim,0.40,0.45,0.55];
    for i = 1:4
        
        myAC.Xcg_adim = X_CG_vector(i);
        myAC.Cm_alpha = -myAC.CL_alpha*(myAC.Xn_adim - myAC.Xcg_adim);

        %  Costanti e condizioni iniziali
        g = 9.81; %Accelerazione di gravità [m/s^2]
        xEG_0 = 0; %[m]
        zEG_0 = -4000; %Altitudine [m]
        V0 = 257.0; %Velocità di volo [m/s];
        q0 = convangvel(0.000,'deg/s','rad/s'); %Velocità angolare di beccheggio
        gamma0 = convang(1.000,'deg','rad'); %Angolo di rampa
        [air_Temp0,sound_speed0,air_pressure0,rho0] = atmosisa(-zEG_0);
 
        %% Processo di minimizzazione della funzione di costo
        %  Valore di tentativo iniziale per il design vector
        x0 = [0;    %Valore di tentativo iniziale per alpha0 [rad]
              0;    %Valore di tentativo iniziale per delta_e0 [rad]
              0;    %Valore di tentativo iniziale per delta_s0 [rad]
              0.5]; %Valore di tentativo iniziale per delta_T0
      
        %  Definizione del vincolo al parametro delta_s0
        Aeq = zeros(4);
        Aeq(3,3) = 1;
        delta_s_0 = convang(-0.5000,'deg','rad');            
        beq = zeros(4,1);
        beq(3,1) = delta_s_0;
        %L'esercizio in esame non richiede di imporre alcun vincolo ai 
        %parametri di trim. Tuttavia, per ottenere i medesimi valori 
        %iniziali deducibili dai diagrammi in figura 7.13, a partire dai 
        %quali realizzare la simulazione del moto, è stato necessario imporre
        %il parametro delta_s0 pari a -0.5deg.
    
        %  Limiti
        lb = [convang(-15,'deg','rad'),... %Valore minimo per alpha
              convang(-20,'deg','rad'),... %Valore minimo per delta_e_0
              convang(-5,'deg','rad'),...  %Valore minimo per delta_s_0
              0.2];                        %Valore minimo per delta_T_0
        ub = [convang(15,'deg','rad'),...  %Valore massimo per alpha
              convang(13,'deg','rad'),...  %Valore massimo per delta_e_0
              convang(5,'deg','rad'),...   %Valore massimo per delta_s_0
              1.0];                        %Valore massimo per delta_T_0

        %  0pzioni di ricerca del minimo
        options = optimset('tolfun',1e-9,'Algorithm','interior-point');

        %  Chiamata alla funzione 'fmincon'
        [x,fval] = fmincon(@costLongEquilibriumStaticStickFixed,...
                           x0,...                                     
                           [],[],Aeq,beq,...                       
                           lb,ub,...                                  
                           @myNonLinearConstraint,...               
                           options);                                  

        alpha0_rad = x(1);
        trim_matrix(1,i) = convang(alpha0_rad,'rad','deg');
        theta0_rad = alpha0_rad - myAC.mu_x + gamma0;
        trim_matrix(2,i) = convang(theta0_rad,'rad','deg');
        delta_e0_rad = x(2);
        trim_matrix(3,i) = convang(delta_e0_rad,'rad','deg');
        delta_s0_rad = x(3);
        trim_matrix(4,i) = convang(delta_s0_rad,'rad','deg');
        delta_T0 = x(4);
        trim_matrix(5,i) = delta_T0;

        %% Integrazione delle equazioni del moto a 3-DoF
        t_fin = 100; %Tempo di simulazione [s]
        state0 = [V0,alpha0_rad,q0,xEG_0,zEG_0,theta0_rad];

        %  Assegnazione delle leggi temporali dei comandi di volo
        global delta_e...
               delta_s...
               delta_T

        delta_e = @(t) interp1([0,t_fin],[delta_e0_rad,delta_e0_rad],t,'linear');
        delta_s = @(t) interp1([0,t_fin],[delta_s0_rad,delta_s0_rad],t,'linear');
        delta_T = @(t) interp1([0,t_fin],[delta_T0,delta_T0],t,'linear');

        %  Integrazione del sistema di equazioni differenziali
        options = odeset('RelTol', 1e-9,'AbsTol', 1e-9*ones(1,6));
        [vTime,mState] = ode45(@eqLongDynamicStickFixed,[0 t_fin],state0,options);
       
        if i == 1
           A = [vTime,mState];
           vDelta_e_A = convang(delta_e(vTime),'rad','deg');
           vDelta_s_A = convang(delta_s(vTime),'rad','deg');
           vDelta_T_A = delta_T(vTime); 
        elseif i == 2
               B = [vTime,mState];
               vDelta_e_B = convang(delta_e(vTime),'rad','deg');
               vDelta_s_B = convang(delta_s(vTime),'rad','deg');
               vDelta_T_B = delta_T(vTime);
        elseif i == 3
               C = [vTime,mState];
               vDelta_e_C = convang(delta_e(vTime),'rad','deg');
               vDelta_s_C = convang(delta_s(vTime),'rad','deg');
               vDelta_T_C = delta_T(vTime);
        elseif i == 4
               D = [vTime,mState];
               vDelta_e_D = convang(delta_e(vTime),'rad','deg');
               vDelta_s_D = convang(delta_s(vTime),'rad','deg');
               vDelta_T_D = delta_T(vTime);
        end

    end
    
end

%% Grafica
%Diagrammi dei comandi di volo al variare di X_CG_adim (figure 1-4)
figure(1)
subplot 311
plot(A(:,1),vDelta_e_A,'b-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([-7 1]) 
ylabel('$\delta_e(deg)$','interpreter','latex','fontsize',11);
subplot 312
plot(A(:,1),vDelta_s_A,'b-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([-1 1]) 
ylabel('$\delta_s(deg)$','interpreter','latex','fontsize',11);
subplot 313
plot(A(:,1),vDelta_T_A,'b-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([0 1])
xlabel('$t (s)$','interpreter','latex','fontsize',11);
ylabel('$\delta_T$','interpreter','latex','fontsize',11);
sgtitle('$\hat{X}_{CG} = 0.32$','interpreter','latex','fontsize',11)

figure(2)
subplot 311
plot(B(:,1),vDelta_e_B,'b-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([-7 1]) 
ylabel('$\delta_e(deg)$','interpreter','latex','fontsize',11);
subplot 312
plot(B(:,1),vDelta_s_B,'b-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([-1 1]) 
ylabel('$\delta_s(deg)$','interpreter','latex','fontsize',11);
subplot 313
plot(B(:,1),vDelta_T_B,'b-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([0 1]) 
xlabel('$t (s)$','interpreter','latex','fontsize',11);
ylabel('$\delta_T$','interpreter','latex','fontsize',11);
sgtitle('$\hat{X}_{CG} = 0.40$','interpreter','latex','fontsize',11)

figure(3)
subplot 311
plot(C(:,1),vDelta_e_C,'b-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([-7 1]) 
ylabel('$\delta_e(deg)$','interpreter','latex','fontsize',11);
subplot 312
plot(C(:,1),vDelta_s_C,'b-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([-1 1]) 
ylabel('$\delta_s(deg)$','interpreter','latex','fontsize',11);
subplot 313
plot(C(:,1),vDelta_T_C,'b-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([0 1]) 
xlabel('$t (s)$','interpreter','latex','fontsize',11);
ylabel('$\delta_T$','interpreter','latex','fontsize',11);
sgtitle('$\hat{X}_{CG} = 0.45$','interpreter','latex','fontsize',11)

figure(4)
subplot 311
plot(D(:,1),vDelta_e_D,'b-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([-7 1]) 
ylabel('$\delta_e(deg)$','interpreter','latex','fontsize',11);
subplot 312
plot(D(:,1),vDelta_s_D,'b-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([-1 1]) 
ylabel('$\delta_s(deg)$','interpreter','latex','fontsize',11);
subplot 313
plot(D(:,1),vDelta_T_D,'b-.','LineWidth',1.5);
grid on
xlim([0 t_fin])
ylim([0 1]) 
xlabel('$t (s)$','interpreter','latex','fontsize',11);
ylabel('$\delta_T$','interpreter','latex','fontsize',11);
sgtitle('$\hat{X}_{CG} = 0.55$','interpreter','latex','fontsize',11)

%Storie temporali delle variabili di stato (figure 5-6-7)
vVel_A = A(:,2);
vAlpha_rad_A = A(:,3);
v_u_A = vVel_A.*cos(vAlpha_rad_A - myAC.mu_x);
v_w_A = vVel_A.*sin(vAlpha_rad_A - myAC.mu_x);
v_q_A = convangvel(A(:,4),'rad/s','deg/s');
vAlpha_deg_A = convang(vAlpha_rad_A,'rad','deg');
vAlphaB_deg_A = vAlpha_deg_A - convang(myAC.mu_x,'rad','deg');
vTheta_deg_A = convang(A(:,7),'rad','deg');
vGamma_deg_A = vTheta_deg_A - vAlphaB_deg_A;
vXe_A = A(:,5);
vZe_A = A(:,6);
vDelta_V_A = vVel_A - V0;
vDelta_alpha_A = vAlpha_deg_A - trim_matrix(1,1);
vDelta_q_A = v_q_A - convangvel(q0,'rad/s','deg/s');
vDelta_h_A = -(vZe_A - zEG_0+V0*sin(gamma0)*A(:,1));

vVel_B = B(:,2);
vAlpha_rad_B = B(:,3);
v_u_B = vVel_B.*cos(vAlpha_rad_B - myAC.mu_x);
v_w_B = vVel_B.*sin(vAlpha_rad_B - myAC.mu_x);
v_q_B = convangvel(B(:,4),'rad/s','deg/s');
vAlpha_deg_B = convang(vAlpha_rad_B,'rad','deg');
vAlphaB_deg_B = vAlpha_deg_B - convang(myAC.mu_x,'rad','deg');
vTheta_deg_B = convang(B(:,7),'rad','deg');
vGamma_deg_B = vTheta_deg_B - vAlphaB_deg_B;
vXe_B = B(:,5);
vZe_B = B(:,6);
vDelta_V_B = vVel_B - V0;
vDelta_alpha_B = vAlpha_deg_B - trim_matrix(1,2);
vDelta_q_B = v_q_B - convangvel(q0,'rad/s','deg/s');
vDelta_h_B = -(vZe_B - zEG_0+V0*sin(gamma0).*B(:,1));

vVel_C = C(:,2);
vAlpha_rad_C = C(:,3);
v_u_C = vVel_C.*cos(vAlpha_rad_C - myAC.mu_x);
v_w_C = vVel_C.*sin(vAlpha_rad_C - myAC.mu_x);
v_q_C = convangvel(C(:,4),'rad/s','deg/s');
vAlpha_deg_C = convang(vAlpha_rad_C,'rad','deg');
vAlphaB_deg_C = vAlpha_deg_C - convang(myAC.mu_x,'rad','deg');
vTheta_deg_C = convang(C(:,7),'rad','deg');
vGamma_deg_C = vTheta_deg_C - vAlphaB_deg_C;
vXe_C = C(:,5);
vZe_C = C(:,6);
vDelta_V_C = vVel_C - V0;
vDelta_alpha_C = vAlpha_deg_C - trim_matrix(1,3);
vDelta_q_C = v_q_C - convangvel(q0,'rad/s','deg/s');
vDelta_h_C = -(vZe_C - zEG_0+V0*sin(gamma0)*C(:,1));

vVel_D = D(:,2);
vAlpha_rad_D = D(:,3);
v_u_D = vVel_D.*cos(vAlpha_rad_D - myAC.mu_x);
v_w_D = vVel_D.*sin(vAlpha_rad_D - myAC.mu_x);
v_q_D = convangvel(D(:,4),'rad/s','deg/s');
vAlpha_deg_D = convang(vAlpha_rad_D,'rad','deg');
vAlphaB_deg_D = vAlpha_deg_D - convang(myAC.mu_x,'rad','deg');
vTheta_deg_D = convang(D(:,7),'rad','deg');
vGamma_deg_D = vTheta_deg_D - vAlphaB_deg_D;
vXe_D = D(:,5);
vZe_D = D(:,6);
vDelta_V_D = vVel_D - V0;
vDelta_alpha_D = vAlpha_deg_D - trim_matrix(1,4);
vDelta_q_D = v_q_D - convangvel(q0,'rad/s','deg/s');
vDelta_h_D = -(vZe_D - zEG_0-V0*sin(gamma0)*D(:,1));

figure(5)
subplot 211
plot(A(:,1),vDelta_V_A,'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5);
hold on;
plot(B(:,1),vDelta_V_B,'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
plot(C(:,1),vDelta_V_C,':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
plot(D(:,1),vDelta_V_D,'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',1.5);
grid on
lgd = legend('$\hat{X}_{CG} = 0.32$','$\hat{X}_{CG} = 0.40$','$\hat{X}_{CG} = 0.45$','$\hat{X}_{CG} = 0.55$');
lgd.Interpreter = 'latex'; 
lgd.FontSize = 11;
xlim([0 t_fin])
ylim([-1e-2 0.15])
ylabel('$\Delta V (m)$','interpreter','latex','fontsize',11)
subplot 212
plot(A(:,1),vDelta_alpha_A,'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5);
hold on;
plot(B(:,1),vDelta_alpha_B,'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
plot(C(:,1),vDelta_alpha_C,':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
plot(D(:,1),vDelta_alpha_D,'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',1.5);
grid on
lgd = legend('$\hat{X}_{CG} = 0.32$','$\hat{X}_{CG} = 0.40$','$\hat{X}_{CG} = 0.45$','$\hat{X}_{CG} = 0.55$');
lgd.Interpreter = 'latex'; 
lgd.FontSize = 11;
xlim([0 t_fin])
xlabel('$t (s)$','interpreter','latex','fontsize',11);
ylim([-0.05e-2 0.5e-2])
ylabel('$\Delta \alpha (deg)$','interpreter','latex','fontsize',11)

figure(6)
subplot 211
plot(A(:,1),vDelta_q_A,'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5);
hold on;
plot(B(:,1),vDelta_q_B,'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
plot(C(:,1),vDelta_q_C,':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
plot(D(:,1),vDelta_q_D,'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',1.5);
grid on
lgd = legend('$\hat{X}_{CG} = 0.32$','$\hat{X}_{CG} = 0.40$','$\hat{X}_{CG} = 0.45$','$\hat{X}_{CG} = 0.55$');
lgd.Interpreter = 'latex'; 
lgd.FontSize = 11;
xlim([0 t_fin])
ylim([-1e-2 1.2e-3])
ylabel('$\Delta q (deg/s)$','interpreter','latex','fontsize',11)
subplot 212
plot(A(:,1),vDelta_h_A,'-.','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5);
hold on;
plot(B(:,1),vDelta_h_B,'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
plot(C(:,1),vDelta_h_C,':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
plot(D(:,1),vDelta_h_D,'-.','color',[0.4940, 0.1840, 0.5560],'LineWidth',1.5);
grid on
lgd = legend('$\hat{X}_{CG} = 0.32$','$\hat{X}_{CG} = 0.40$','$\hat{X}_{CG} = 0.45$','$\hat{X}_{CG} = 0.55$');
lgd.Interpreter = 'latex'; 
lgd.FontSize = 11;
xlim([0 t_fin])
ylim([-1 0.001])
xlabel('$t (s)$','interpreter','latex','fontsize',11);
ylabel('$\Delta h (m)$','interpreter','latex','fontsize',11)

figure(7)
subplot 221
plot(A(:,1),vTheta_deg_A,'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
hold on;
plot(A(:,1),vGamma_deg_A,':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
grid on
lgd = legend('$\theta(t)$','$\gamma(t)$');
lgd.Interpreter = 'latex'; 
lgd.FontSize = 11;
xlim([0 t_fin])
%ylim([-0.01 0.08])
ylabel('$(deg)$','interpreter','latex','fontsize',11)
title('$\hat{X}_{CG} = 0.32$','interpreter','latex','fontsize',11)
subplot 222
plot(B(:,1),vTheta_deg_B,'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
hold on;
plot(B(:,1),vGamma_deg_B,':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
grid on
lgd = legend('$\theta(t)$','$\gamma(t)$');
lgd.Interpreter = 'latex'; 
lgd.FontSize = 11;
xlim([0 t_fin])
%ylim([-4e-2 1e-2])
ylabel('$(deg)$','interpreter','latex','fontsize',11)
title('$\hat{X}_{CG} = 0.40$','interpreter','latex','fontsize',11)
subplot 223
plot(C(:,1),vTheta_deg_C,'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
hold on;
plot(C(:,1),vGamma_deg_C,':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
grid on
lgd = legend('$\theta(t)$','$\gamma(t)$');
lgd.Interpreter = 'latex'; 
lgd.FontSize = 11;
xlim([0 t_fin])
xlabel('$t (s)$','interpreter','latex','fontsize',11);
%ylim([-10e-2 2e-2])
ylabel('$(deg)$','interpreter','latex','fontsize',11)
title('$\hat{X}_{CG} = 0.45$','interpreter','latex','fontsize',11)
subplot 224
plot(D(:,1),vTheta_deg_D,'--','color',[0, 0.4470, 0.7410],'LineWidth',1.5);
hold on;
plot(D(:,1),vGamma_deg_D,':','color',[0.4660, 0.6740, 0.1880],'LineWidth',1.5);
grid on
lgd = legend('$\theta(t)$','$\gamma(t)$');
lgd.Interpreter = 'latex'; 
lgd.FontSize = 11;
xlim([0 t_fin])
xlabel('$t (s)$','interpreter','latex','fontsize',11);
ylim([-5 65])
ylabel('$(deg)$','interpreter','latex','fontsize',11)
title('$\hat{X}_{CG} = 0.55$','interpreter','latex','fontsize',11)