clear all; clc; close all;

global myAC g rho0 q0 gamma0 z0 deltas0 deltae0 deltaT0 delta_e delta_T delta_s

%% Input properties
aircraftDataFileName = 'DSV_Aircraft_data.txt';
myAC = DSVAircraft(aircraftDataFileName);
g = 9.81;

%% Initial Trim Conditions
aircraftDataFileName = 'DSV_Aircraft_data.txt';
myAC = DSVAircraft(aircraftDataFileName);
tf = 300;

if (myAC.err == -1)
    disp('Termination.')
else
    % Constants and initial conditions
    xEG_0 = 0; % [m]
    zEG_0 = -4000; % Altitude [m]
    q0 = convangvel(0.000, 'deg/s', 'rad/s'); % Pitch angular velocity
    gamma0 = convang(0.00, 'deg', 'rad'); % Ramp angle
    [air_Temp0, sound_speed0, air_pressure0, rho0] = atmosisa(-zEG_0);
    v0 = 0.5 * sound_speed0; %Mach 0.5

    % Process of cost function minimization
    x0 = [0; 0; 0; 0.5]; % Initial guess for the design vector
    Aeq = zeros(4);
    Aeq(3, 3) = 1;
    delta_s_0 = convang(-1.000, 'deg', 'rad');
    beq = zeros(4, 1);
    beq(3, 1) = delta_s_0;%ho fissato delta_s_0

    lb = [convang(-15, 'deg', 'rad'), convang(-20, 'deg', 'rad'), convang(-5, 'deg', 'rad'), 0.2];
    ub = [convang(15, 'deg', 'rad'), convang(13, 'deg', 'rad'), convang(5, 'deg', 'rad'), 1.0];

    options = optimset('tolfun', 1e-9, 'Algorithm', 'interior-point');

    global V_0 q_0 gamma_0 rho_0
    V_0 = v0;
    q_0 = q0;
    gamma_0 = gamma0;
    rho_0 = rho0;

    [x, fval] = fmincon(@(x) costLongEquilibriumStaticStickFixed(x), x0, ...
        [], [], Aeq, beq, lb, ub, @myNonLinearConstraint, options);

    alpha0 = x(1);
    alpha_0_deg = convang(x(1), 'rad', 'deg');
    delta_e_0 = x(2);
    deltae0_deg = convang(x(2), 'rad', 'deg');
    deltas0 = x(3);
    delta_s_0_deg = convang(x(3), 'rad', 'deg');
    deltaT0 = x(4);
    theta0 = alpha0 + gamma0 - myAC.mu_x;
    z0 = -4000;
    xE0 = 0;
    [~, ~, ~, rho0] = atmosisa(-z0);

    vBreaksdelta_e(1, :) = [0, 4.99, 5, 6.0, 6.01, 14.99, 15, 16.0, 16.01, tf];
    vBreaksdelta_e(2, :) = [deltae0_deg, deltae0_deg, deltae0_deg-1,...
        deltae0_deg-1, deltae0_deg, deltae0_deg, deltae0_deg+1, deltae0_deg+1, ...
        deltae0_deg, deltae0_deg];

    delta_e_deg = @(t) interp1(vBreaksdelta_e(1, :), vBreaksdelta_e(2, :), t, 'linear');
    delta_e = @(t) convang(delta_e_deg(t), 'deg', 'rad');
    delta_T = @(t) interp1 ([0 tf], [deltaT0 deltaT0], t);
    delta_s = @(t) interp1 ([0 tf], [deltas0 deltas0], t);
    t_span=linspace(0,tf,1000); %vettore dei tempi
    x0 = [v0, alpha0, q0, xE0, z0, theta0];%condizioni iniziali
    options = odeset('RelTol', 1e-3, 'AbsTol', 1e-3 * ones(1, 6));
    [vTime, mState] = ode45(@eqLongDynamicStickFixed, t_span, x0);
    vdeltae = delta_e_deg(vTime);
end

%%%%%%%%%Parte linearizzata%%%%%%%%%%
h = -zEG_0;
Mach_0 = .5; %Mach
alpha_b = alpha_0_deg;
% h_0 = convlength(h,'ft','m');
U_0 = Mach_0 * sound_speed0*cos(convang(alpha_b,'deg','rad')); %è uguale alla V_0 di prima
%
S = myAC.S; % m^2
cbar = myAC.mac; % m
Weight_0 = myAC.W; % N
mass = Weight_0/g;
Iyy_0 = mass*myAC.k_y^2;
Iyy = Iyy_0;

C_L = 2*Weight_0/(rho0*U_0^2*S);
C_D = myAC.CD_0 + myAC.K*(C_L)^myAC.m;

C_L_alpha = myAC.CL_alpha; % 1/rad
C_D_alpha = 2*myAC.CL_alpha*alpha0; %Tramite la polare
C_m_alpha_0 = myAC.Cm_alpha; % 1/rad
C_L_de = myAC.CL_delta_e; % 1/rad
C_m_de = myAC.Cm_delta_e; % 1/rad
C_L_alphadot = myAC.CL_alpha_dot; % 1/rad
C_m_alphadot = myAC.Cm_alpha_dot; % 1/rad
C_L_q = myAC.CL_q; % 1/rad
C_m_q = myAC.Cm_q; % 1/rad
C_L_Mach = 0 ; %trascurato
C_D_Mach = 0; %trascurato
C_m_Mach = 0.10; %valore plausibile

%% Data elaboration
% KV = 0; CTFIX = 0;
%% Constants Computation -----> Quaderno 16, pg. 60
qbar_0 = 0.5*rho0*(U_0^2);
mu_0 = mass/(0.5*rho0*S*cbar);
Gamma_0 = 0.0; % rad
SM_0 = 0.22;
SM = 0.22;
C_m_alpha = C_m_alpha_0*(SM/SM_0); % 1/rad
X_u = -(qbar_0*S/(mass*U_0))*(2*C_D + Mach_0*C_D_Mach); % constant thrust
X_w = (qbar_0*S/(mass*U_0))*(C_L - C_D_alpha);
X_wdot = 0; X_q = 0;
X_de = 0;
X_dT = qbar_0*S/mass*(0 + 0/U_0^2); %
Z_u = -(qbar_0*S/(mass*U_0))*( ...
    2*C_L + (Mach_0^2/(1-Mach_0^2))*C_L_Mach); % constant thrust
Z_w = -(qbar_0*S/(mass*U_0))*(C_D + C_L_alpha); Z_wdot = -(1/(2*mu_0))*C_L_alphadot;
Z_q = -(U_0/(2*mu_0))*C_L_q; Z_de = -qbar_0*S/mass*C_L_de;
Z_dT = 0;
%
M_u = (qbar_0*S*cbar/(Iyy*U_0))*Mach_0*C_m_Mach; M_w = (qbar_0*S*cbar/(Iyy*U_0))*C_m_alpha;
M_wdot = (rho0*S*(cbar^2)/(4*Iyy))*C_m_alphadot; M_q = (rho0*U_0*S*(cbar^2)/(4*Iyy))*C_m_q;
M_de = qbar_0*S/Iyy*C_m_de; M_dT = 0;
%
% Z_de =Z_de*100; M_de = M_de*100;
k_hat = M_wdot/(1-Z_wdot);
% Plant matrix --> cfr. (16.147b)
A_lon(1,1) = X_u;
A_lon(1,2) = X_w; A_lon(1,3) = 0;
A_lon(1,4) = -g*cos(Gamma_0); A_lon(2,1) = Z_u/(1 - Z_wdot);
A_lon(2,2) = Z_w/(1 - Z_wdot); A_lon(2,3) = (Z_q + U_0)/(1 - Z_wdot);
A_lon(2,4) = -g*sin(Gamma_0)/(1 - Z_wdot); A_lon(3,1) = M_u + k_hat*Z_u;
A_lon(3,2) = M_w + k_hat*Z_w; A_lon(3,3) = M_q + k_hat*(Z_q+U_0);
A_lon(3,4) = -k_hat*g*sin(Gamma_0); A_lon(4,1) = 0;
A_lon(4,2) = 0; A_lon(4,3) = 1;
A_lon(4,4) = 0;
%
%chiamo eig che prende in ingresso la matrice quadrata e resittuisce una
%matrice di autovalori V e autovalori D posi sulla diagonale
[V,D] = eig(A_lon);
%
W = inv(V); %ha gli autovalori sulle righe
% X_de/mass, X_dT/mass; ...
% Z_de/(mass-Z_wdot), Z_dT/(mass-Z_wdot); ...
% (M_de + Z_de*M_wdot/(mass-Z_wdot))/Iyy, ...
% (M_dT + Z_dT*M_wdot/(mass-Z_wdot))/Iyy; ...
% 0, 0];
B_lon = [ X_de , X_dT; ...
    Z_de/(1-Z_w) , Z_dT/(1-Z_w); ...
    M_de + k_hat* Z_de , M_dT + k_hat* Z_dT ; ...
    0, 0];
% Output matrices
C_lon = eye(3,4);
D_lon = zeros(3,2);
% make a state-space representation
sys_lon = ss( ...
    A_lon , ... % A
    B_lon , ... % B
    eye(4,4), ... % C
    zeros(4,2) ... % D
    );
% % Time vector
t_BP = [0, 4.99, 5, 6.0, 6.01, 14.99, 15, 16.0, 16.01,tf];
de_deg_BP = [0, 0 , -1, -1 , 0 , 0 , 1, 1, 0 , 0];

% figure; plot(t_BP, de_deg_BP); % Time vector
t_Elevator_Doublet = linspace(0,tf,1000)';

% Commanded deflection column vector

de_rad_Elevator_Doublet = interp1(t_BP,convang(de_deg_BP,'deg','rad'),t_Elevator_Doublet);
de_deg_Elevator_Doublet = interp1(t_BP,de_deg_BP,t_Elevator_Doublet); % Input u, nx2 matrix

u_Elevator_Doublet = [ ...
    de_rad_Elevator_Doublet, zeros(length(de_rad_Elevator_Doublet),1)];

% Initial condition
x0_Elevator_Doublet = [0;0;0;0];
% Simulation
[y_Elevator_Doublet , t_Elevator_Doublet , x_Elevator_Doublet] = ...
    lsim(sys_lon , u_Elevator_Doublet , t_Elevator_Doublet , x0_Elevator_Doublet);
%La funzione lsim restituisce la risposta del sistema y, campionata negli 
% stessi istanti temporali t dell'input. Per sistemi con un solo output, 
% y è un vettore della stessa lunghezza di t. Per sistemi con più output,
% y è una matrice con tante righe quante sono le misurazioni nel tempo 
% (lunghezza di t) e tante colonne quante sono le uscite nel sistema sys.
% Questa sintassi non genera un grafico.
% x_elevator_doublet contiene le perturbazioni per colonne ovvero 1000 sono
% le riche corrispondono al numero dei tempoi e lungo le colonne ci sono le
% variazioni di alpha, w , q e theta 
W_0 = U_0*tand(alpha_b);
Vel_0 = sqrt(U_0^2 + W_0^2);

d_V = sqrt((U_0+x_Elevator_Doublet(:,1)).^2+ ...
    (W_0+x_Elevator_Doublet(:,2)).^2)- sqrt(U_0.^2+W_0.^2);

d_alpha = atan((W_0+x_Elevator_Doublet(:,2))./(U_0+x_Elevator_Doublet(:,1)))-atan( W_0/U_0);

d_alpha_deg = convang(d_alpha,'rad','deg'); %

d_q_rad_s = x_Elevator_Doublet(:,3); d_q_deg_s = convangvel(d_q_rad_s,'rad/s','deg/s' );

d_theta_rad = x_Elevator_Doublet(:,4);

% figure; plot(t_Elevator_Doublet , alpha_0_deg+d_alpha_deg);
%% Comparative plots
figure; plot(vTime,vdeltae, 'b-',LineWidth=1.5); grid on; hold on;
plot(t_Elevator_Doublet,deltae0_deg + de_deg_Elevator_Doublet, 'r--',LineWidth=1.5 );
title({'Deflessione dell''equilibratore'},'fontsize',15);legend('Q7','Q16'); ylabel('\delta_e'); xlabel('time (s)');
figure; subplot(2,1,1);
plot(vTime,convang(mState(:,2),'rad','deg')); grid on; hold on;
plot(t_Elevator_Doublet , alpha_0_deg+d_alpha_deg);title({'Angolo d''attacco'},'fontsize',15);
xlabel('t (s)'); ylabel('\alpha (deg)'); legend('Q7', 'Q16');
subplot(2,1,2);
plot(vTime,mState(:,1)); grid on; hold on;
plot(t_Elevator_Doublet , Vel_0+d_V);
title({'Velocità'},'fontsize',15);
xlabel('t (s)'); ylabel('V (m/s)');
legend('Q7', 'Q16');
figure; subplot(2,1,1);
plot(vTime,convangvel(mState(:,3),'rad/s','deg/s')); grid on; hold on;
plot(t_Elevator_Doublet , d_q_deg_s );title({'Velocità angolare di beccheggio'},'fontsize',15);
xlabel('t (s)'); ylabel('q (deg/s)'); legend('Q7', 'Q16')
subplot(2,1,2);
plot(vTime,convang(mState(:,6),'rad','deg')); grid on; hold on;
plot(t_Elevator_Doublet, convang(theta0 + d_theta_rad,'rad','deg'));
title({'Angolo di elevazione'},'fontsize',15);
xlabel('t (s)'); ylabel('\theta (deg)'); legend('Q7', 'Q16')
