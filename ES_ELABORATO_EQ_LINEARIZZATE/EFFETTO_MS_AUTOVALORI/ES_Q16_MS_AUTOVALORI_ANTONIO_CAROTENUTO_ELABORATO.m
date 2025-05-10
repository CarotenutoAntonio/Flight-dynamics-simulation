%Definizioni delle grandezze da inserire nelle matrici
%%%% metto qua poi correggo%%%%%%
clc; clear; clear all;
% Grandezze geometriche
condition = 3;
if condition==3
SM_vec=[0.22 , 0.15 0.10 , 0.05 , 0.035 , 0.025 , 0.0158 , 0.008 , ...
    0.005 , 0.0021 , 0.0 , -0.002 , -0.005 , -0.008 , -0.0145];
lvec=length(SM_vec);
elseif condition==4
    SM_vec=[0.22 , 0.15 0.10 , 0.05 , 0.035 , 0.025 , 0.0258]% , 0.008 , ...
    %0.005 , 0.0021 , 0.0 , -0.002 , -0.005 , -0.008 , -0.0145];
    lvec=length(SM_vec);
    elseif condition==1
    SM_vec=[0.22 , 0.15 0.10 , 0.05 , 0.035 , 0.025 , 0.0258]% , 0.008 , ...
    %0.005 , 0.0021 , 0.0 , -0.002 , -0.005 , -0.008 , -0.0145];
    lvec=length(SM_vec);
    elseif condition==2
    SM_vec=[0.22 , 0.15 0.10 , 0.05 , 0.035 , 0.025 , 0.0258]% , 0.008 , ...
    %0.005 , 0.0021 , 0.0 , -0.002 , -0.005 , -0.008 , -0.0145];
    lvec=length(SM_vec);
elseif condition==5
    SM_vec=[0.22 , 0.15 0.10 , 0.05 , 0.035 , 0.025 , 0.0258]% , 0.008 , ...
    %0.005 , 0.0021 , 0.0 , -0.002 , -0.005 , -0.008 , -0.0145];
    lvec=length(SM_vec);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ac data
mass = [255753 2.8869e+05 2.8869e+05 2.8869e+05 2.8869e+05]; %vettore delle masse
mass= mass( condition ); %kg
Iyy = [4.38*10^7 4.2740*10^7 4.2740*10^7 4.2740*10^7 4.2740*10^7];%vettore momenti di inerzia
Iyy = Iyy ( condition ) ;
S = 510.97;%superificie alare
cbar = 8.32;% corda media aerodinamica
SM_0 = 0.22;%margine statico di sicurezza iniziale
%% Condizioni di volo
zEG_0 = [0 2e+4 2e+4 4e+4 4e+4]; %quote di volo
zEG_0 = zEG_0*0.3048;% Conversione ft—>m
zEG_0 = zEG_0( condition) ;
q0 = 0;% Velocità angolare di beccheggio
Gamma_0 = 0;% Angolo di rampa
g_0 = 9.81; %g
[~, a_0 ,~, rho_0] = atmosisa ( zEG_0 ) ; %ISA
Mach_0 = [0.25 0.5 0.8 0.8 0.9];% vettore Mach
Mach_0 = Mach_0( condition ) ;
U_0 = Mach_0* a_0 ; %velocità
alfa_B_0 = [5.70 6.80 0.00 4.60 2.40] ;%vettore alpha body
alfa_B_0 = alfa_B_0( condition ) ;
qbar_0 = 0.5* rho_0*U_0^2; %pressione dinamica
mu_0 = 2*mass/(rho_0 *S*cbar ) ;

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
C_m_alpha_0 = [-1.26 -1.15 -0.63 -1.03 -1.61];
C_m_alpha_0=C_m_alpha_0(condition);


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



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Derivate di stabilità
 
     X_u = -(qbar_0*S/(mass*U_0))*(2*C_D + Mach_0*C_D_Mach);                 % Dipende da D e T quindi Cd ed M
    X_w = (qbar_0*S/(mass*U_0))*(C_L - C_D_alpha);                         % dipende da W->L->Cl
    Z_u = -(qbar_0*S/(mass*U_0))*(2*C_L + (Mach_0^2/(1-Mach_0^2))*C_L_Mach); % constant thrust
    Z_w = -(qbar_0*S/(mass*U_0))*(C_D + C_L_alpha);                        % Proporzionale a -CLalfa
    Z_wdot = -(1/(2*mu_0))*C_L_alphadot;
    Z_q = -(U_0/(2*mu_0))*C_L_q;                                           
    M_u = (qbar_0*S*cbar/(Iyy*U_0))*Mach_0*C_m_Mach;
                           
    M_wdot = (rho_0*S*(cbar^2)/(4*Iyy))*C_m_alphadot;                       % Proviene dal downwash
    M_q = (rho_0*U_0*S*(cbar^2)/(4*Iyy))*C_m_q;                             % Proporzionale a CMq (derivata di smorzamento)(-)
    k_hat = M_wdot/(1-Z_wdot);      



    %%%%%%%%%% messo qua il ciclo
    for i=1:lvec
    SM=SM_vec(i);                 
    C_m_alpha=C_m_alpha_0*(SM/SM_0);

    M_w    = (qbar_0*S*cbar/(Iyy*U_0))... % [1/(ms)]
        *C_m_alpha;


    A_lon(1,1) = X_u;
    A_lon(1,2) = X_w;
    A_lon(1,3) = 0;
    A_lon(1,4) = -g_0*cos(Gamma_0);              % Prima riga

    A_lon(2,1) = Z_u/(1 - Z_wdot);
    A_lon(2,2) = Z_w/(1 - Z_wdot);
    A_lon(2,3) = (Z_q + U_0)/(1 - Z_wdot);
    A_lon(2,4) = -g_0*sin(Gamma_0)/(1 - Z_wdot); % Seconda riga

    A_lon(3,1) = M_u + k_hat*Z_u;
    A_lon(3,2) = M_w + k_hat*Z_w;
    A_lon(3,3) = M_q + k_hat*(Z_q+U_0);
    A_lon(3,4) = -k_hat*g_0*sin(Gamma_0);        % Terza riga

    A_lon(4,1) = 0;
    A_lon(4,2) = 0;
    A_lon(4,3) = 1;
    A_lon(4,4) = 0;                              % Quarta riga
    [V,D] = eig(A_lon);
    W = inv(V);
    eigen_vals_vec{i,1}= D ;
    eigen_vals_vec{i,2}= D ;
end
for i=1:lvec
    lambda_SP_arr(i)=eigen_vals_vec{i,1}(1 ,1);
    lambda_SP2_arr(i)=eigen_vals_vec{i,1}(2,2);
    sigma_SP_arr(i)=real(lambda_SP_arr(i));
    omega_SP_arr(i)=imag(lambda_SP_arr(i));
    sigma_SP2_arr(i)=real(lambda_SP2_arr(i));
    omega_SP2_arr(i)=imag(lambda_SP2_arr(i));

    lambda_PH_arr(i)=eigen_vals_vec{i,1}(3,3);
    lambda_PH2_arr(i)=eigen_vals_vec{i,1}(4,4);
    sigma_PH_arr(i)=real(lambda_PH_arr(i));
    omega_PH_arr(i)=imag(lambda_PH_arr(i));
    sigma_PH2_arr(i)=real(lambda_PH2_arr(i));
    omega_PH2_arr(i)=imag(lambda_PH2_arr(i));
end

%%Plots
figure(1);
plot(sigma_PH_arr, omega_PH_arr, 'ro','LineWidth', 1, 'MarkerSize', 8); hold on;
plot(sigma_PH2_arr, omega_PH2_arr, 'bo','LineWidth', 1, 'MarkerSize', 8);
plot(sigma_SP_arr, omega_SP_arr, 'go','LineWidth', 1, 'MarkerSize', 8);
plot(sigma_SP2_arr, omega_SP2_arr, 'ko','LineWidth', 1, 'MarkerSize', 8);
xline(0, '-k', 'LineWidth', 0.5);  % Linea dell'asse x in nero tratteggiato
yline(0, '-k', 'LineWidth', 0.5);  % Linea dell'asse y in nero tratteggiato
legend('Phugoid 1', 'Phugoid 2','Short Period 1', 'Short Period 2');
title('Eigenvalues for different SMs');
xlabel( 'Re');
ylabel('Im');

