%% ESERCIZIO 16_1: DA UNO SCRIPT MATLAB RICAVARE LA MATRICE (16.196) CALCOLANDONE AUTOVALORI ED AUTOVETTORI...
% DIAGRAMMA RISPOSTE MODALI E FATTORI DEGLI AUTOVETTORI. VERIFICA ESISTENZA MODI CORTO PERIODO E FUGOIDE E CALCOLARNE...
% LE CARATTERISTICHE: COEFF SMORZ, PULS NAT, PERIODO, TEMPO DI DIMEZZ, NUMERO CICLI FINO AL TEMPO DI DIMEZZ.

clc; clear all; close all;
%%
    %% modifico qua in modo da poter considerare diverse condizioni
%%

%% Scelta della condizione di volo
condition = 3;
%% ac data
mass = [255753 2.8869e+05 2.8869e+05 2.8869e+05 2.8869e+05]; %vettore delle masse
mass= mass( condition ); %kg
Iyy = [4.38*10^7 4.2740*10^7 4.2740*10^7 4.2740*10^7 4.2740*10^7];%vettore momenti di inerzia
Iyy = Iyy ( condition ) ;
S = 510.97;%superificie alare
cbar = 8.32;% corda media aerodinamica
SM0 = 0.22;%margine statico di sicurezza iniziale
SM= 0.22 ;%(X_N - X_G) / cbar

%% Condizioni di volo
zEG_0 = [0 2e+4 2e+4 4e+4 4e+4]; %quote di volo
zEG_0 = zEG_0*0.3048;% Conversione ft—>m
zEG_0 = zEG_0( condition) ;
q0 = 0;% Velocità angolare di beccheggio
gamma0 = 0;% Angolo di rampa
g = 9.81; %g
[~, a_0 ,~, rho0] = atmosisa ( zEG_0 ) ; %ISA
Mach0 = [0.25 0.5 0.8 0.8 0.9];% vettore Mach
Mach0 = Mach0( condition ) ;
U_0 = Mach0* a_0 ; %velocità
alfa_B_0 = [5.70 6.80 0.00 4.60 2.40] ;%vettore alpha body
alfa_B_0 = alfa_B_0( condition ) ;
qbar_0 = 0.5* rho0*U_0^2; %pressione dinamica
mu_0 = 2*mass/(rho0 *S*cbar ) ;

%Definizione dei vettori delle caratteristiche aerodinamiche e delle derivate
% di stabilità
C_L = [1.10 0.68 0.27 0.66 0.52];
C_L =C_L (condition);
C_D = [0.10 0.04 0.02 0.04 0.04];
C_D = C_D( condition ) ;
C_L_alpha = [5.70 4.67 4.24 4.92 5.57];
C_L_alpha =C_L_alpha (condition);
C_D_alpha = [0.66 0.37 0.08 0.43 0.53];
C_D_alpha=C_D_alpha(condition);
c_m_alpha_0 = [-1.26 -1.15 -0.63 -1.03 -1.61];
c_m_alpha_0=c_m_alpha_0(condition);
C_m_alpha = c_m_alpha_0 *SM/SM0;

C_L_alphadot = [6.70 6.53 5.99 5.91 5.53];
C_L_alphadot = C_L_alphadot( condition) ;
C_m_alphadot = [-3.20 3.35 -5.40 -6.41 -8.82];
C_m_alphadot= C_m_alphadot (condition ) ;

C_L_q = [5.40 5.13 5.01 6.00 6.94];
C_L_q = C_L_q( condition ) ;
C_m_q = [-20.80 -20.70 -20.50 -24 -25.10];
C_m_q = C_m_q( condition ) ;
C_L_Mach = [0 -0.09 0.11 0.21 -0.28] ;
C_L_Mach = C_L_Mach( condition ) ;
C_D_Mach = [0 0 0.01 0.03 0.24];
C_D_Mach = C_D_Mach( condition ) ;
C_m_Mach = [0 0.12 -0.12 0.17 -0.11] ;
C_m_Mach = C_m_Mach( condition ) ;
C_L_de = [0.338 0.356 0.270 0.367 0.300];
C_L_de = C_L_de ( condition ) ;
C_m_de = [-1.34 -1.43 1.06 -1.45 -1.20];
C_m_de = C_m_de ( condition ) ;

% Calcolo derivare di stabilità 
% Deriv stab longitudinali
    X_u = -(qbar_0*S/(mass*U_0))*(2*C_D + Mach0*C_D_Mach);                 % Dipende da D e T quindi Cd ed M
    X_w = (qbar_0*S/(mass*U_0))*(C_L - C_D_alpha);                         % dipende da W->L->Cl
    Z_u = -(qbar_0*S/(mass*U_0))*(2*C_L + (Mach0^2/(1-Mach0^2))*C_L_Mach); % constant thrust
    Z_w = -(qbar_0*S/(mass*U_0))*(C_D + C_L_alpha);                        % Proporzionale a -CLalfa
    Z_wdot = -(1/(2*mu_0))*C_L_alphadot;
    Z_q = -(U_0/(2*mu_0))*C_L_q;                                           
    M_u = (qbar_0*S*cbar/(Iyy*U_0))*Mach0*C_m_Mach;
    M_w = (qbar_0*S*cbar/(Iyy*U_0))*C_m_alpha;                             % dipende da W->L->Cm_alpha
    M_wdot = (rho0*S*(cbar^2)/(4*Iyy))*C_m_alphadot;                       % Proviene dal downwash
    M_q = (rho0*U_0*S*(cbar^2)/(4*Iyy))*C_m_q;                             % Proporzionale a CMq (derivata di smorzamento)(-)
    k_hat = M_wdot/(1-Z_wdot);      




% Costruzione matrice  A_LON
    A = NaN(4);
    A(1,1) = X_u;
    A(1,2) = X_w;
    A(1,3) = 0;
    A(1,4) = -g*cos(gamma0);
    A(2,1) = Z_u/(1-Z_wdot);
    A(2,2) = Z_w/(1-Z_wdot);
    A(2,3) = (Z_q + U_0)/(1-Z_wdot);
    A(2,4) = -g*sin(gamma0)/(1-Z_wdot);
    A(3,1) = M_u + k_hat*Z_u;
    A(3,2) = M_w + k_hat*Z_w;
    A(3,3) = M_q + k_hat*(Z_q + U_0);
    A(3,4) = -k_hat*g*sin(gamma0);
    A(4,1) = 0;
    A(4,2) = 0;
    A(4,3) = 1;
    A(4,4) = 0;

    % Derivate di controllo
    X_dT = 0;
    Z_dT = 0;
    M_dT = 0;
    X_de = 0;
    Z_de = -qbar_0*S/mass*C_L_de;
    M_de = qbar_0*S*cbar/Iyy*C_m_de;

% RIPOSTA LIBERA E CARATTERISTICHE MODALI 
% determinazione autovalori e autovettori
% Troviamo autovalori ed autovettori del sistema: V ci restituisce una
% matrice di autovalori e D la matrice dei corrispondenti autovettori tale
% che A_lon * V = V * D
    [V,D] = eig(A);

% Per avere una scala dei valori si dividono tutte le componenti per la
% quarta componente tale che V_SP(4) = V_fugoide(4) = 1 + 0*i.
% I valori della seconda e quarta colonna non vengono considerati solo
% perché sono semplicemente i complessi coniugati ripettivamente della
% prima e terza colonna

% Autovalori corto periodo rispetto a theta (4)
    V_SP = V(:,1);
    V_SP = V_SP/V_SP(4);
    V_SP(1) = V_SP(1)/U_0;
    V_SP(2) = V_SP(2)/U_0;
    V_SP(3) = V_SP(3)/(2*U_0/cbar);

% Autovalori lungo periodo rispetto a theta (4)
    V_Ph = V(:,3);
    V_Ph = V_Ph/V_Ph(4);
    V_Ph(1) = V_Ph(1)/U_0;
    V_Ph(2) = V_Ph(2)/U_0;
    V_Ph(3) = V_Ph(3)/(2*U_0/cbar);


% SS ANALIZZA SPAZIO DEGLI STATI (PER SISTEMA LINEARE TEMPO INVARIANTE)
% Definizione matrice
    sys = ss( ...
    A, ...           % A
    zeros(4,1), ...  % B
    eye(4,4), ...    % C
    zeros(4,1) ...   % D
    );

% CORTO PERIODO E FUGOIDE 
    V1SP = V(:,1);
    V1SP  = V1SP/V1SP(4);
    V1PH = V(:,3);
    V1PH  = V1PH/V1PH(4);

% Valori iniziali
    X0_1 = real(V1SP);
    X0_2 = real(V1PH);

% Risuluzione problema LTI ai valori iniziali
% (INITIAL) VALUTA SOLUZIONI PROBLEMA (AUTOVALORI)
    [Y_1,T_1,X_1] = initial(sys,X0_1);
    [Y_2,T_2,X_2] = initial(sys,X0_2);


% ESERCIZIO 16.2
% Calcolo caratteristiche modali richieste
    disp('CORTO PERIODO');
    zita_SP     = sqrt(1/(1+(imag(D(1,1))/(real(D(1,1))))^2));
    disp(['smorzamento ',num2str(zita_SP)]);
    omega_n_SP = -real(D(1,1))/(zita_SP);
    disp(['frequenza ',num2str(omega_n_SP),' s^-1']);
    T_SP        = (2*pi)/(omega_n_SP*sqrt(1-zita_SP^2));
    disp(['durata ',num2str(T_SP),' s']);
    t_half_SP   = (log(2))/(omega_n_SP*zita_SP);
    disp(['tempo di dimezzamento ',num2str(t_half_SP),'s']);
    N_half_SP   = t_half_SP/T_SP;
    disp(['numero di cicli di dimezzamento ',num2str(N_half_SP),' s']);

disp(' ');

disp('FUGOIDE');
zita_Ph     = sqrt(1/(1+(imag(D(3,3))/(real(D(3,3))))^2));
disp(['smorzamento ',num2str(zita_Ph)]);
omega_n_Ph = -real(D(3,3))/(zita_Ph);
disp(['frequenza ',num2str(omega_n_Ph),' s^-1']);
T_Ph        = (2*pi)/(omega_n_Ph*sqrt(1-zita_Ph^2));
disp(['durata ',num2str(T_Ph),' s']);
t_half_Ph   = (log(2))/(omega_n_Ph*zita_Ph);
disp(['tempo di dimezzamento ',num2str(t_half_Ph),' s']);
N_half_Ph   = t_half_Ph/T_Ph;
disp(['numero di cicli di dimezzamento ',num2str(N_half_Ph),' s']);

disp(' ');
disp(' ');

% PLOT CORTO PERIODO
figure(1);
c = compass(V_SP);
c1 = c(1);
c1.LineWidth = 2;
c1.Color = 'r';

c2 = c(2);
c2.LineWidth = 2;
c2.Color = 'g';


c3 = c(3);
c2.LineWidth = 2;
c3.Color = 'b';

c4 = c(4);
c4.LineWidth = 2;
c4.Color = 'k';


axis('square');
title({'CORTO PERIODO'},'fontsize',15);
legend('$\frac{u}{U_0}$', '$\frac{w}{U_0}$', '$\frac{q c}{2U_0}$', '$\theta$',...
    'FontSize',15);
set(legend,...
    'Location','bestoutside',...
    'Orientation','vertical',...
    'Interpreter','latex',...
    'FontSize',15);



% PLOT FUGOIDE
figure(2);
cc = compass(V_Ph);
c1 = cc(1);
c1.LineWidth = 2;
c1.Color = 'r';

c2 = cc(2);
c2.LineWidth = 2;
c2.Color = 'g';


c3 = cc(3);
c2.LineWidth = 2;
c3.Color = 'b';

c4 = cc(4);
c4.LineWidth = 2;
c4.Color = 'k';


axis('square');
title({'FUGOIDE'},'fontsize',15);
legend('$\frac{u}{U_0}$', '$\frac{w}{U_0}$', '$\frac{q c}{2U_0}$', '$\theta$',...
    'FontSize',15);
set(legend,...
    'Location','bestoutside',...
    'Orientation','vertical',...
    'Interpreter','latex',...
    'FontSize',15);



figure(3); 
plot(T_1, Y_1(:,1)/U_0, 'r','linewidth',1.2); % u/U0
hold on;
plot(T_1, Y_1(:,2)/U_0, 'g','linewidth',1.2); % w/U0
hold on;
plot(T_1, Y_1(:,3)*((cbar)/(2*U_0)), 'b','linewidth',1.2); % qc/2U0
hold on 
plot(T_1,Y_1(:,4), 'k','linewidth',1.2); % theta
title({'CORTO PERIODO'},'fontsize',15);
legend('$\frac{u}{U_0}$', '$\frac{w}{U_0}$', '$\frac{q c}{2U_0}$', '$\theta$',...
    'FontSize',15);
set(legend,...
    'Location','best',...
    'Orientation','horizontal',...
    'Interpreter','latex',...
    'FontSize',15);

xlabel('t (s)', 'FontSize',12);
ylim([-1 1])
xlim([0 T_1(end)])
grid on;


figure(4);
plot(T_2, Y_2(:,1)/U_0, 'r','linewidth',1.2);% u/U0
hold on;
plot(T_2, Y_2(:,2)/U_0, 'g','linewidth',1.2);% w/U0
hold on;
plot(T_2, Y_2(:,3)*((cbar)/(2*U_0)),'b','linewidth',1.2);% qc/2U0
hold on 
plot(T_2,Y_2(:,4),'k','linewidth',1.2);% theta
title({'FUGOIDE'},'fontsize',15);
legend('$\frac{u}{U_0}$', '$\frac{w}{U_0}$', '$\frac{q c}{2U_0}$', '$\theta$',...
    'FontSize',15);
set(legend,...
    'Location','best',...
    'Orientation','horizontal',...
    'Interpreter','latex',...
    'FontSize',15);
xlim([0 T_2(end)/3])
xlabel('t (s)', 'FontSize',12);
grid on;

% GRAFICA AUTOVALORI PIANO COMPLESSO
figure(5); 
plot(real(D(1,1)),imag(D(1,1)),'k.','markersize',20);
hold on;
plot(real(D(2,2)),imag(D(2,2)),'k.','markersize',20);
hold on;
plot(real(D(3,3)),imag(D(3,3)),'b.','markersize',20);
hold on;
plot(real(D(4,4)),imag(D(4,4)),'b.','markersize',20);
title({'RAPPRESENTAZIONE AUTOVALORI PIANO Re-Im'},'fontsize',15);
grid on;
xlim([-2.5 0.5])
xlabel('Re','FontSize',12)
ylabel('Im','FontSize',12)
legend('\lambda_{SP}', '\lambda_{SP}^*', '\lambda_{PH}',...
    '\lambda_{PH}^*',...
    'Location','best',...
    'fontsize',12);



% ESERCIZIO 16.6
% CARATTERISTICHE APPROSSIMATE DI CORTO PERIODO E FUGOIDE
% APPROSSIMAZIONE DI CORTO PERIODO
    disp('CARATTERISTICHE APPROSSIMATE DI CORTO PERIODO');
    omega_n_SP = sqrt(Z_w*M_q - M_w*U_0);
    disp(['frequenza ',num2str(omega_n_SP),' s^-1']);
    zita_SP    = -(Z_w +M_q + M_wdot*U_0)/(2*omega_n_SP);
    disp(['smorzamento ',num2str(zita_SP)]);
    T_SP       = (2*pi)/(omega_n_SP*sqrt(1-zita_SP^2));
    disp(['durata ',num2str(T_SP),' s']);
    t_half_SP  = (log(2))/(omega_n_SP*zita_SP);
    disp(['tempo di dimezzamento ',num2str(t_half_SP),' s']);
    N_half_SP  = t_half_SP/T_SP;
    disp(['numero di cicli di dimezzamento ',num2str(N_half_SP),' s']);

disp(' ');

% APPROSSIMAZIONE DI CORTO PERIODO "coarse"
    disp('CARATTERISTICHE APPROSSIMATE DI CORTO PERIODO "coarse"');
    omega_n_SPc = sqrt(- M_w*U_0);
    disp(['frequenza ',num2str(omega_n_SPc),' s^-1']);
    zita_SPc    = -(M_q)/(2*omega_n_SPc);
    disp(['smorzamento ',num2str(zita_SPc)]);
    T_SPc       = (2*pi)/(omega_n_SPc*sqrt(1-zita_SPc^2));
    disp(['durata ',num2str(T_SPc),' s']);
    t_half_SPc  = (log(2))/(omega_n_SPc*zita_SPc);
    disp(['tempo di dimezzamento ',num2str(t_half_SPc),' s']);
    N_half_SPc  = t_half_SPc/T_SP;
    disp(['numero di cicli di dimezzamento ',num2str(N_half_SPc),' s']);

disp(' ');

% APPROSSIMAZIONE DI FUGOIDE
A_approx = [ X_u + X_w*((M_u*U_0 - M_q*Z_u)/(-M_w*U_0 + M_q*Z_q)),  -g;...
                             (M_w*Z_u-M_u*Z_w)/(-M_w*U_0+M_q*Z_w),   0];
[V_approx,D_approx]      = eig(A_approx);    
    disp('CARATTERISTICHE APPROSSIMATE DI FUGOIDE');
    zita_PH    = sqrt(1/(1+(imag(D_approx(1,1))/(real(D_approx(1,1))))^2));
    disp(['smorzamento ',num2str(zita_PH)]);
    omega_n_PH = -real(D_approx(1,1))/(zita_PH);
    disp(['frequenza ',num2str(omega_n_PH),' s^-1']);
    T_PH       = (2*pi)/(omega_n_PH*sqrt(1-zita_PH^2));
    disp(['durata ',num2str(T_PH),' s']);
    t_half_PH  = (log(2))/(omega_n_SP*zita_PH);
    disp(['tempo di dimezzamento  ',num2str(t_half_PH),'s']);
    N_half_PH  = t_half_PH/T_PH;
    disp(['numero di cicli di dimezzamento ',num2str(N_half_PH),' s']);

    % APPROSSIMAZIONE DI FUGOIDE "coarse"
A_approxc = [ X_u ,  -g;...
                             (-Z_u)/U_0,   0];
[V_approxc,D_approxc]      = eig(A_approxc);    
    disp('CARATTERISTICHE APPROSSIMATE DI FUGOIDE "coarse"');
    zita_PHc    = sqrt(1/(1+(imag(D_approxc(1,1))/(real(D_approxc(1,1))))^2));
    disp(['smorzamento ',num2str(zita_PHc)]);
    omega_n_PHc = -real(D_approxc(1,1))/(zita_PHc);
    disp(['frequenza ',num2str(omega_n_PHc),' s^-1']);
    T_PHc       = (2*pi)/(omega_n_PHc*sqrt(1-zita_PHc^2));
    disp(['durata ',num2str(T_PHc),' s']);
    t_half_PHc  = (log(2))/(omega_n_SP*zita_PHc);
    disp(['tempo di dimezzamento  ',num2str(t_half_PHc),'s']);
    N_half_PHc  = t_half_PHc/T_PHc;
    disp(['numero di cicli di dimezzamento ',num2str(N_half_PHc),' s']);
    





