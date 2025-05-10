clc;
clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CONDIZIONE 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_fin = 100;
% Inerzie 
Weight = convforce(636640,'lbf','N');
g_0    = 9.81;
mass   = Weight/g_0;
S       = 5500*(convlength(1,'ft','m'))^2;% Superficie alare
h=0;
[T_0, a_0, P_0, rho_0]=atmosisa(h);
cbar    = convlength(27.3,'ft','m');       % mac                 [m]
mu_0    = mass/(0.5*rho_0*S*cbar);
Gamma_0 = 0.0;                             % angolo di salita iniziale [rad]
AR=9.45;
b=sqrt(S*AR);
M=0.25;
U_0=a_0*M;
qbar_0  = 0.5*rho_0*(U_0^2);               % pressione dinamica iniziale [N/m^2]
SLUGFT2toKGM2 = convmass(1,'slug','kg')...   % Si passa da [slug*ft^2] a
    *(convlength(1,'ft','m')^2); % [kg*m^2]
Ixx = 14.30e+6*SLUGFT2toKGM2 ; %14.30e+6 * 0.3048^2; % kg*m^2
Izz = 45.30e+6*SLUGFT2toKGM2 ;% kg*m^2
Ixz = -2.23e+6 *SLUGFT2toKGM2 ;% kg*m^2
i1 = Ixz / Ixx; %prodotti di inerzia
i2 = Ixz / Izz;

% Derivate di stabilità Latero-Direzionale 
Clbeta = -0.221;
Clp = -0.450;
Clr = 0.101;
Cybeta = -0.96;
Cyp = 0;
Cyr = 0.0;
Cnbeta = 0.150;
CnTbeta = 0;
Cnp = -0.121;
Cnr = -0.300;

% derivate di controllo Latero-Direzionale 
Cldeltaa = 0.0461;
Cldeltar = 0.007;
Cydeltaa = 0;
Cydeltar = 0.175;
Cndeltaa = 0.0064;
Cndeltar = -0.109;

% derivate di stabilità dimensionale Latero-Direzionale
Ybeta = (qbar_0 * S * Cybeta) / mass;
Yp = (qbar_0 * S * b * Cyp) / (2 * mass * U_0);
Yr = (qbar_0 * S * b * Cyr) / (2 * mass * U_0);
Lbeta = (qbar_0 * S * b * Clbeta) / Ixx;
Lp = (qbar_0 * S * b^2 * Clp) / (2 * Ixx * U_0);
Lr = (qbar_0 * S * b^2 * Clr) / (2 * Ixx * U_0);
Nbeta = (qbar_0 * S * b * Cnbeta) / Izz;
NTbeta = (qbar_0 * S * b * CnTbeta) / Izz;
Np = (qbar_0 * S * b^2 * Cnp) / (2 * Izz * U_0);
Nr = (qbar_0 * S * b^2 * Cnr) / (2 * Izz * U_0);

%derivate di controllo  Latero-Direzionalel
Ydeltaa = (qbar_0 * S * Cydeltaa) / mass;
Ydeltar = (qbar_0 * S * Cydeltar) / mass;
Ldeltaa = (qbar_0 * S * b * Cldeltaa) / Ixx;
Ldeltar = (qbar_0 * S * b * Cldeltar) / Ixx;
Ndeltaa = (qbar_0 * S * b * Cndeltaa) / Izz;
Ndeltar = (qbar_0 * S * b * Cndeltar) / Izz;

% derivate prime di stabilità Latero-Direzionale
Ybeta_1 = Ybeta;
Yp_1 = Yp;
Yr_1 = Yr;
Lbeta_1 = (Lbeta + i1 * Nbeta) / (1 - i1 * i2);
Lp_1 = (Lp + i1 * Np) / (1 - i1 * i2);
Lr_1 = (Lr + i1 * Nr) / (1 - i1 * i2);
Nbeta_1 = (i2 * Lbeta + Nbeta) / (1 - i1 * i2);
Np_1 = (i2 * Lp + Np) / (1 - i1 * i2);
Nr_1 = (i2 * Lr + Nr) / (1 - i1 * i2);

% Derivate prime di controllo Latero-Direzionali 
Ydeltaa_1 = Ydeltaa;
Ydeltar_1 = Ydeltar;
Ldeltaa_1 = (Ldeltaa + i1 * Ndeltaa) / (1 - i1 * i2);
Ldeltar_1 = (Ldeltar + i1 * Ndeltar) / (1 - i1 * i2);
Ndeltaa_1 = (i2 * Ldeltaa + Ndeltaa) / (1 - i1 * i2);
Ndeltar_1 = (i2 * Ldeltar + Ndeltar) / (1 - i1 * i2);

% Matrice Latero-Direzionale 
A_ld = [Nr_1, Nbeta_1, Np_1, 0;...
    Yr_1/U_0 - 1, Ybeta_1/U_0, Yp_1/U_0, g_0/U_0;...
    Lr_1, Lbeta_1, Lp_1, 0;...
    0, 0, 1, 0];

B_ld = [Ndeltaa_1, Ndeltar_1;...
    Ydeltaa_1/U_0, Ydeltar_1/U_0;...
    Ldeltaa_1, Ldeltar_1;...
    0, 0];

%comando eig prende come input la matrice quadrata e dà come output la
%matrice Vld che contine gli autovettori e Dld gli autovalori sulla
%diagonale 
[Vld, Dld] = eig(A_ld);
%estraggo gli autovettori
V_ld_DR = Vld(:, 2);
% V_ld_DR = V_ld_DR / V_ld_DR(4, 1);
V_ld_Spiral = Vld(:, 4);
% V_ld_Spiral = V_ld_Spiral / V_ld_Spiral(4, 1);
V_ld_Roll = Vld(:, 1);
% V_ld_Roll = V_ld_Roll / V_ld_Roll(4, 1);
sys = ss(...
    A_ld, ... % A
    B_ld, ... % B
    eye(4,4), ... % C
    zeros(4,2) ... % D
);
x0 = [0; 0; 0; 0];
time_dense = [0:0.25:t_fin]';
u_null = [0 * time_dense, 0 * time_dense];
%si poteva utilizzare anche la funzione initial come nel caso longitudinale

[y_DR, t_DR, x_DR] = lsim(sys, u_null, time_dense, real(V_ld_DR));
[y_Roll, t_Roll, x_Roll] = lsim(sys, u_null, time_dense, real(V_ld_Roll));
[y_Spiral, t_Spiral, x_Spiral] = lsim(sys, u_null, time_dense, real(V_ld_Spiral));

% Fasori
V_ld_DR = V_ld_DR / V_ld_DR(4, 1);
V_ld_Spiral = V_ld_Spiral / V_ld_Spiral(4, 1);
V_ld_Roll = V_ld_Roll / V_ld_Roll(4, 1);

% Colori
my_blue = [0, 0.4470, 0.7410];
my_orange = [0.8500, 0.3250, 0.0980];
my_yellow = [0.9290, 0.6940, 0.1250];
my_violet = [0.4940, 0.1840, 0.5560];
my_green = [0.4660, 0.6740, 0.1880];
my_cyan = [0.3010, 0.7450, 0.9330];
my_red = [0.6350, 0.0780, 0.1840];

figure 
title('Dutch Roll Dynamics');
grid on;
plot(t_DR,x_DR(:,1),'b',LineWidth=1);
xlim([0 60]);
hold on;
plot(t_DR,x_DR(:,2),'r',LineWidth=1);
hold on;
xlim([0 60]);
plot(t_DR,x_DR(:,3),'y',LineWidth=1);
hold on;
xlim([0 60]);
plot(t_DR,x_DR(:,4),'m',LineWidth=1);
hold on;
xlim([0 60]);
legend('$r$', '$\beta$', '$p$', '$\phi$', 'Interpreter', 'Latex', 'Location', 'northeastoutside');
title('Dutch Roll Dynamics');
figure
grid on;
plot(t_DR,x_Roll(:,1),'b',LineWidth=1);
xlim([0 10]);
hold on;
plot(t_DR,x_Roll(:,2),'r',LineWidth=1);
hold on;
xlim([0 10]);
plot(t_DR,x_Roll(:,3),'y',LineWidth=1);
hold on;
xlim([0 10]);
plot(t_DR,x_Roll(:,4),'m',LineWidth=1);
hold on;
xlim([0 10]);
legend('$r$', '$\beta$', '$p$', '$\phi$', 'Interpreter', 'Latex', 'Location', 'northeastoutside');
title(' Roll Dynamics');
figure 
grid on;
plot(t_DR,x_Spiral(:,1),'b',LineWidth=1);
xlim([0 100]);
hold on;
plot(t_DR,x_Spiral(:,2),'r',LineWidth=1);
hold on;
xlim([0 100]);
plot(t_DR,x_Spiral(:,3),'y',LineWidth=1);
hold on;
xlim([0 100]);
plot(t_DR,x_Spiral(:,4),'m',LineWidth=1);
hold on;
xlim([0 100]);
legend('$r$', '$\beta$', '$p$', '$\phi$', 'Interpreter', 'Latex', 'Location', 'northeastoutside');
title('Spiral Dynamics');

% Figure
fig_Q16_1_Phasor_Roll = figure(4);
V_roll_3 = compass(V_ld_Roll(3)); hold on;
V_roll_4 = compass(V_ld_Roll(4)); hold on;
V_roll_2 = compass(V_ld_Roll(2)); hold on;
V_roll_1 = compass(V_ld_Roll(1)); hold on;

% Plot
W.Color = [1 1 1];
plotROLL = [V_roll_1 V_roll_2 V_roll_3 V_roll_4];
V_roll_1.Color = my_cyan;
V_roll_1.LineWidth = 1;
V_roll_2.Color = my_blue; 
V_roll_2.LineWidth = 2; % r
V_roll_3.Color = my_yellow;
V_roll_3.LineWidth = 1; % PHI
V_roll_4.Color = my_orange; 
V_roll_4.LineWidth = 1;

lgd = legend(plotROLL, '$r$', '$\beta$', '$p$', '$\phi$', 'Interpreter', ...
    'Latex', 'Location', 'northeastoutside');
lgd.FontSize = 20;
title('Roll Phasors');
set(gcf, 'position', [278.6, 0.2, 600, 400]); % aspect ratio 3/2

% Figure per fasori
fig_Q16_B737_Phasor_DR = figure(5);

% fasori per DR
V_DR_3 = compass(V_ld_DR(3)); hold on; % p
V_DR_4 = compass(V_ld_DR(4)); hold on; % PHI
V_DR_2 = compass(V_ld_DR(2)); hold on; % BETA
V_DR_1 = compass(V_ld_DR(1)); hold on; % r

% Plot
plotPH_DR = [V_DR_1 V_DR_2 V_DR_3 V_DR_4];
V_DR_1.Color = my_cyan; V_DR_1.LineWidth = 2;
V_DR_2.Color = my_blue; V_DR_2.LineWidth = 2;
V_DR_3.Color = my_yellow; V_DR_3.LineWidth = 2;
V_DR_4.Color = my_orange; V_DR_4.LineWidth = 2;

lgd_DR = legend(plotPH_DR, '$r$', '$\beta$', '$p$', '$\phi$', 'Interpreter', 'Latex', 'Location', 'northeastoutside');
lgd_DR.FontSize = 20;
title('Dutch-Roll Phasors');
set(gcf, 'position', [278.6, 0.2, 600, 400]); % aspect ratio 3/2

% Figure per  fasori di spirale
fig_Q16_B737_Phasor_Spiral = figure(6);

% Fasori per Spirale
V_Spiral_4 = compass(V_ld_Spiral(4)); hold on;
V_Spiral_2 = compass(V_ld_Spiral(2)); hold on;
V_Spiral_3 = compass(V_ld_Spiral(3)); hold on;
V_Spiral_1 = compass(V_ld_Spiral(1)); hold on;

% Plot
plotPH_Spiral = [V_Spiral_1 V_Spiral_2 V_Spiral_3 V_Spiral_4];
V_Spiral_1.Color = my_cyan; V_Spiral_1.LineWidth = 2;
V_Spiral_2.Color = my_blue; V_Spiral_2.LineWidth = 2;
V_Spiral_3.Color = my_yellow; V_Spiral_3.LineWidth = 2;
V_Spiral_4.Color = my_orange; V_Spiral_4.LineWidth = 2;

lgd_Spiral = legend(plotPH_Spiral, '$r$', '$\beta$', '$p$', '$\phi$', ...
    'Interpreter', 'Latex', 'Location', 'northeastoutside');
lgd_Spiral.FontSize = 20;
title('Spiral Phasors');
set(gcf, 'position', [278.6, 0.2, 600, 400]); % aspect ratio 3/2




