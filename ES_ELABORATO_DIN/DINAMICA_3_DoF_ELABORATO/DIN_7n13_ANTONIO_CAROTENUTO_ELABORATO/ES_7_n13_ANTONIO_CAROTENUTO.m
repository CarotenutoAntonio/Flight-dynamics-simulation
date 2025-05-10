clc;
clear;
close all;

global...
Adata ... % Aircraft data. Velivolo completo 
WBflap ... % Wing-body con flap 
v_alpha ... % vettore angolo di attacco 
v_delta_e ... % vettore deflessione equilibratore 
v_delta_flap...% vettore deflessione flaps 
v_Mach ... % vettore numeri di Mach 
Alpha_M ... % vettore angoli alpha al variare con il Mach 
Mach_a ... % vettore numero di mach al variare di a 
Deltae_M ... % vettore deflessione equilibratore con il Mach 
Mach_de ... % Mach al variare della deflessione dell'equi. 
Deltaf_M ... % vettore deflessione flap con il numero di Mach 
Mach_df ... % vettore mach al variare della deflessione flap 
Alpha_de_M ... % 
Deltae_a_M ... % 
Mach_a_de ... % 
Alpha_df_M ... % 
Deltaf_a_M ... % 
Mach_a_df ... %


alldata=datcomimport('CITATION_ES_7n13_RIPROVO.out');

WBflap=alldata{1,1};
Adata=alldata{1,2};


%Vettori delle variabili di input 
v_alpha = convang(Adata.alpha, 'deg', 'rad'); % [rad] 
v_delta_e = convang(Adata.delta, 'deg', 'rad'); % [rad] 
v_delta_flap = convang(WBflap.delta, 'deg', 'rad'); % [rad] 
v_Mach = Adata.mach;

%Creazione delle matrici dei nodi di supporto 

[Mach_a, Alpha_M] = meshgrid(v_Mach ,v_alpha); %% per Cma,Cla,CL,CM,CD(clear)
[Mach_de, Deltae_M] = meshgrid(v_Mach ,v_delta_e); %% per Dcl_de,Dcm_de,Dcdmin_de
[Mach_df, Deltaf_M] = meshgrid(v_Mach ,v_delta_flap); %% per Dcl_df,Dcm_df,Dcdmin_df
[Alpha_de_M, Deltae_a_M, Mach_a_de] = ndgrid(v_alpha ,v_delta_e ,v_Mach);%% per Dcdi_de
[Alpha_df_M, Deltaf_a_M, Mach_a_df] = ndgrid(v_alpha ,v_delta_flap ,v_Mach); %% per Dcdi_df

% provo le funzioni
%[a,b]=calcola_Cmadot_Cladot_ANTONIO_CAROTENUTO(10.5/57.3,0.5);
%[c,d,e]=calcola_CL_CM_CD_ANTONIO_CAROTENUTO(10.5/57.3,-10/57.3,25/57.3,0.5);
%% FINE COSTRUZIONE COEFFICIENTI




%% INIZIO PROBLEMA DI DINAMICA

disp('Moto del velivolo a 3 gradi di libertà');
disp('Risoluzione del problema di trim ad una data altitudine e velocità di volo');

%% Dichiarazione delle variabili globali
global g...                   %Accelerazione di gravità 
       zEG_0 V0 q0 gamma0...  %Condizioni iniziali
       rho0 ...               %Densità dell'aria all'altitudine h = (-zEG_0)
       myAC                   %Oggetto 'Velivolo'

%% Ricerca delle condizioni di trim
%aircraftDataFileName = 'DSV_Aircraft_data_ES7_n13.txt';

%  Definizione dell'oggetto 'Velivolo'
%myAC = DSVAircraft(aircraftDataFileName);
%L'oggetto 'myAC' è un'istanza della classe DSVAircraft.
%Tutte le variabili membro di 'myAC' risultano essere accessibili 
%e modificabili dall'utente.
%%
%%%ho modificato qua
%%%
%%

%% Ricerca delle condizioni di trim
aircraftDataFileName = 'DSV_Aircraft_data.txt';

%  Definizione dell'oggetto 'Velivolo'
myAC = DSVAircraft(aircraftDataFileName);


if (myAC.err == -1)
    disp('Terminazione.')
else
    disp(['File ',aircraftDataFileName,' letto correttamente.']);   

    % Costanti e condizioni iniziali
    g = 9.81; %Accelerazione di gravità [m/s^2]
    xEG_0 = 0; %[m]
    zEG_0 = -4000; %Altitudine [m]
    V0 = 257.0; %Velocità di volo [m/s]
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
    lb = [convang(-12,'deg','rad'),... %Valore minimo per alpha
          convang(-20,'deg','rad'),... %Valore minimo per delta_e0
          convang(-5,'deg','rad'),...  %Valore minimo per delta_s0
          0.2];                        %Valore minimo per delta_T0
    ub = [convang(11,'deg','rad'),...  %Valore massimo per alpha
          convang(13,'deg','rad'),...  %Valore massimo per delta_e0
          convang(5,'deg','rad'),...   %Valore massimo per delta_s0
          1.0];                        %Valore massimo per delta_T0


        %  Definizione del vincolo imposto al parametro delta_s_0
    Aeq = zeros(4);
    Aeq(3,3) = 1;
    delta_s_0 = convang(1.000,'deg','rad');            
    beq = zeros(4,1);
    beq(3,1) = delta_s_0;


    %  0pzioni di ricerca del minimo
    options = optimset('tolfun',1e-9,'Algorithm','interior-point');

    %  Chiamata alla funzione 'fmincon'
    [x,fval] = fmincon(@costLESSF_DATCOM_ANTONIO_CAROTENUTO,... %Funzione obiettivo del processo di minimizzazione
                       x0,...                                     
                       [],[],Aeq,beq,... %Vincoli lineari
                       lb,ub,...                                  
                       @myNonLinearConstraint,... %Vincoli non lineari
                       options);                                  
    %La funzione 'myNonLinearConstraint' fornisce due matrici atte alla
    %definizione di vincoli non lineari. La sua definizione è resa 
    %necessaria anche in assenza di tali vincoli. 
    
    %  Stampa delle componenti del design vector
    alpha0_rad=x(1);
    delta_e0_rad=x(2);
    delta_s0_rad=x(3);

    alpha0_deg = convang(x(1),'rad','deg');
    delta_e0_deg = convang(x(2),'rad','deg');
    delta_s0_deg = convang(x(3),'rad','deg');
    delta_T0 = x(4);
    
    fprintf('alpha_0 = %f deg. \n',alpha0_deg);
    fprintf('delta_e_0 = %f deg. \n',delta_e0_deg);
    fprintf('delta_s_0 = %f deg. \n',delta_s0_deg);
    fprintf('delta_T_0 = %f deg. \n',delta_T0);
    
    %% inizia esercizio 7.8 note le condizioni di trim
    % devo integrare le equazioni della dinamica per moto 3DOF
    
    % definisco valori comandati dopo l'istante iniziale

    global delta_e... %equilibratore
        delta_flap...   %flap
        delta_T...    %manetta
    
    tfin=100;
%% rispetto all'esercizio 8 modifico qua
for i=1:2
    if i==1
    t1=1;
    t2=2.5;
    t3=4;
    
    elseif i==2
    t1=1;
    t2=1.5;
    t3=2;
    end
    var_delta_e=convang(3,'deg','rad');
    delta_e_comm=delta_e0_rad-var_delta_e;



    delta_e = @(t) interp1([0,t1,t2,t3,tfin],[delta_e0_rad,delta_e0_rad,delta_e_comm,delta_e0_rad,delta_e0_rad],t,'linear');
    delta_flap = @(t) interp1([0,tfin],[delta_s0_rad,delta_s0_rad],t,'linear');
    delta_T = @(t) interp1([0,tfin],[delta_T0,delta_T0],t,'linear');
   
    V00=V0;
    q00=q0;
    zEG_00=zEG_0;
    theta0=gamma0+alpha0_rad-myAC.mu_x;
    state0=[V00, alpha0_rad , q00 , xEG_0 , zEG_00 , theta0];
    options = odeset( ...
    'RelTol', 1e-9, ...
    'AbsTol', 1e-9*ones(1,6) ...
    );
    [vTime,mstate] = ode45(@eqLDSF_Datcom,[0 tfin],state0,options);
    %[vTime,mstate] = ode45(@eqLongDynamicStickFixed_ANTONIO_CAROTENUTO,[0 tfin],state0,options);
    mstate_dot = zeros(length(vTime),6);

for j = 1:length(vTime)
    
    mstate_dot(j,:) = eqLDSF_Datcom(vTime(j),mstate(j,:));
    
end

%% calcolo dei fattori di carico
vgamma= mstate(:,6)-mstate(:,2)+myAC.mu_x;
vgamma_dot=mstate_dot(:,6)-mstate_dot(:,2);
vV= mstate(:,1);
vV_dot= mstate_dot(:,1);
f_xa=-(sin(vgamma)+vV_dot/g);
f_za=(cos(vgamma)+vV.*vgamma_dot/g);

    %grafica dell'evoluzione dinamica dello stato
    figure(1)
    subplot 311
    plot(vTime,mstate(:,1)-V00,'-o','LineWidth',1.5,'MarkerSize',1);
    hold on
    ylabel('$\Delta V(m/s)$','Interpreter','latex');

    grid on

        subplot 312
    plot(vTime,convang(mstate(:,2)-alpha0_rad,'rad','deg'),'-o','LineWidth',1.5,'MarkerSize',1);
    hold on
    ylabel('$\Delta \alpha(deg)$','Interpreter','latex');
    grid on

        subplot 313
    plot(vTime,convangvel(mstate(:,3)-q00,'rad/s','deg/s'),'-o','LineWidth',1.5,'MarkerSize',1);
    hold on
    ylabel('$\Delta q$','Interpreter','latex');
    grid on

        figure(2)
    subplot 311
    plot(vTime,mstate(:,4),'-o','LineWidth',1.5,'MarkerSize',1);
    hold on
    ylabel('$x_{E,G}$','Interpreter','latex');
    grid on

        subplot 312
    plot(vTime,mstate(:,5)-zEG_00,'-o','LineWidth',1.5,'MarkerSize',1);
    hold on
    ylabel('$\Delta z_{E,G} (m)$','Interpreter','latex');
    grid on

        subplot 313
    plot(vTime,convang(mstate(:,6)-theta0,'rad','deg'),'-o','LineWidth',1.5,'MarkerSize',1);
    hold on
    ylabel('$\Delta \theta $','Interpreter','latex');
    grid on

    figure(3)
    subplot 311
    plot([0,t1,t2,t3,tfin],[delta_e0_rad,delta_e0_rad,delta_e_comm,delta_e0_rad,delta_e0_rad]);
    hold on
    ylabel('$\delta_{e}$','Interpreter','latex');
    grid on

    subplot 312
    plot([0,tfin],[delta_s0_rad,delta_s0_rad]);
    hold on
    ylabel('$ \delta_{e}$','Interpreter','latex');  
    grid on

    subplot 313
    plot([0,tfin],[delta_T0,delta_T0]);
    hold on
    ylabel('$\delta_{e}$','Interpreter','latex');

    grid on


    figure(4)
    subplot 211
    plot(vTime,f_xa,'-o','LineWidth',1.5,'MarkerSize',1);
    hold on
    ylabel('$f_{xa}$','Interpreter','latex');
    grid on

    subplot 212
    plot(vTime,f_za,'-o','LineWidth',1.5,'MarkerSize',1);
    hold on
    ylabel('$f_{za}$','Interpreter','latex');
    grid on
end
end