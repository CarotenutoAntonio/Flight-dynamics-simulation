%% Ricerca di condizioni di trim per un moto a 3-DoF
clear all; close all; clc;

disp('Moto del velivolo a 3 gradi di libertà');
disp('Risoluzione del problema di trim ad una data altitudine e velocità di volo');

%% Dichiarazione delle variabili globali
global g...                   %Accelerazione di gravità 
       zEG_0 V0 q0 gamma0...  %Condizioni iniziali
       rho0 ...               %Densità dell'aria all'altitudine h = (-zEG_0)
       myAC                   %Oggetto 'Velivolo'

%% Ricerca delle condizioni di trim
aircraftDataFileName = 'DSV_Aircraft_data.txt';

%  Definizione dell'oggetto 'Velivolo'
myAC = DSVAircraft(aircraftDataFileName);
%L'oggetto 'myAC' è un'istanza della classe DSVAircraft.
%Tutte le variabili membro di 'myAC' risultano essere accessibili 
%e modificabili dall'utente.

if (myAC.err == -1)
    disp('Terminazione.')
else
    disp(['File ',aircraftDataFileName,' letto correttamente.']);   

    % Costanti e condizioni iniziali
    g = 9.81; %Accelerazione di gravità [m/s^2]
    xEG_0 = 0; %[m]
    zEG_0 = -4000; %Altitudine [m]
    Nvel=20;
    vV0 = linspace(220.0,350.0,Nvel); %Velocità di volo [m/s]
    q0 = convangvel(0.000,'deg/s','rad/s'); %Velocità angolare di beccheggio
    gamma0 = convang(0.000,'deg','rad'); %Angolo di rampa
    [air_Temp0,sound_speed0,air_pressure0,rho0] = atmosisa(-zEG_0);

    %% Processo di minimizzazione della funzione di costo
    %  Valore di tentativo iniziale per il design vector
    x0 = [0;    %Valore di tentativo iniziale per alpha_0 espresso in rad
          0;    %Valore di tentativo iniziale per delta_e_0 espresso in rad
          0;    %Valore di tentativo iniziale per delta_s_0 espresso in rad
          0.5]; %Valore di tentativo iniziale per delta_T_0
    %Il valore di tentativo iniziale viene dedotto considerando una condizione
    %media tra i limiti inferiore e superiore.
    
    %  Limiti
    lb = [convang(-15,'deg','rad'),... %Valore minimo per alpha
          convang(-20,'deg','rad'),... %Valore minimo per delta_e0
          convang(-5,'deg','rad'),...  %Valore minimo per delta_s0
          0.2];                        %Valore minimo per delta_T0
    ub = [convang(15,'deg','rad'),...  %Valore massimo per alpha
          convang(13,'deg','rad'),...  %Valore massimo per delta_e0
          convang(5,'deg','rad'),...   %Valore massimo per delta_s0
          1.0];                        %Valore massimo per delta_T0

    %  0pzioni di ricerca del minimo
    options = optimset('tolfun',1e-9,'Algorithm','interior-point');

    %  Chiamata alla funzione 'fmincon'
valpha0_deg=zeros(Nvel,1);
vdelta_e0_deg=zeros(Nvel,1);
vdelta_T0=zeros(Nvel,1);
vdelta_s0_deg=zeros(Nvel,1);

% constrain lineari
delta_s0_rad=convang(1,'deg','rad');
Aeq=zeros(4);
Aeq(3,3)=1;
Beq=zeros(4,1);
Beq(3,1)=delta_s0_rad;

for i=1:Nvel

    V0=vV0(i);
    [x,fval] = fmincon(@costLongEquilibriumStaticStickFixed,... %Funzione obiettivo del processo di minimizzazione
                       x0,...                                     
                       [],[],Aeq,Beq,... %Vincoli lineari
                       lb,ub,...                                  
                       @myNonLinearConstraint,... %Vincoli non lineari
                       options);                                  
    %La funzione 'myNonLinearConstraint' fornisce due matrici atte alla
    %definizione di vincoli non lineari. La sua definizione è resa 
    %necessaria anche in assenza di tali vincoli. 
    
    %  Stampa delle componenti del design vector                      
    valpha0_deg(i) = convang(x(1),'rad','deg');
    vdelta_e0_deg(i) = convang(x(2),'rad','deg');
    vdelta_s0_deg(i) = convang(x(3),'rad','deg');
    vdelta_T0(i) = x(4);

end

%% grafica

figure(1)
subplot(3,1,1);
plot(vV0,valpha0_deg','b-o','LineWidth',1.5);
title('alfa');
ylabel('$\alpha_{0} (deg)$','Interpreter','latex');
xlabel('$V_{0}$','Interpreter','latex');

subplot(3,1,2);
plot(vV0,vdelta_e0_deg','b-o','LineWidth',1.5);
title('\delta_{e,0} (deg)');
ylabel('$\delta_{e,0} (deg)$','Interpreter','latex');
xlabel('$V_{0}$','Interpreter','latex');


subplot(3,1,3);
plot(vV0,vdelta_T0','b-o','LineWidth',1.5);
title('\delta_{T,0} (deg)');
ylabel('$\delta_{T,0} (deg)$','Interpreter','latex');
xlabel('$V_{0}$','Interpreter','latex');

end
 %Results.vVel0 = vV0;
    %Results.alpha0_deg = valpha0_deg;
    %Results.delta_e0_deg = vdelta_e0_deg;
    %Results.delta_T0 = vdelta_T0;
    
    %save('linear_resultsANTONIOCAROTENUTO.mat','Results')