function f = eqLongDynamicStickFixed(t,x)
% Funzione per l'implementazione delle equazioni 7.61

% Dichiarazione variabili globali
global ...
     g ...                    % accelerazione di gravità
     delta_e ...              % legge di comando elevatore
     delta_s ...              % legge di comando stabilizzatore
     delta_T ...              % legge di comando manetta
     myAC                     % oggetto aereo rimpito fuori da questa funzione


% Si assegna un nome delle componenti del vettore di design
V       = x(1);
alpha   = x(2);
q       = x(3);
x_EG    = x(4);
z_EG    = x(5);
theta   = x(6);

% Funzione matlab del modello ISA. Lo usiamo per il calcolo della densità
% Serve per aggiornare i parametri atmosferici dato che z_EG può variare
[air_Temp_0, sound_speed_0, air_pressure_0, rho] = atmosisa(-z_EG);

% densità relativa aereo
mu_rel = (myAC.W/g)/(rho*myAC.S*myAC.b); % m/(rho*S*b)

% Secondi membri delle equazioni dinamiche per volo stazionario
% Si vedano le eq. (7.61)
% eq. (7.61a)
F1 = g.*((delta_T(t)*myAC.T/myAC.W)*cos(alpha - myAC.mu_x + myAC.mu_T) ...
         -sin(myAC.mu_x + theta - alpha) ...
         -((rho*V^2)/(2*(myAC.W/myAC.S))) ...
             *(myAC.CD_0 + myAC.K*((myAC.CL_alpha*alpha ...
                    + myAC.CL_delta_e*delta_e(t) ...
                    + myAC.CL_delta_s*delta_s(t))^myAC.m)));

% eq. (7.61b)
F2 = (( 1. - myAC.CL_q*(myAC.mac/myAC.b)/(4*mu_rel))*q ...
       -(g/V)*(delta_T(t)*myAC.T/myAC.W)...
       *sin(alpha - myAC.mu_x + myAC.mu_T) ...
       +(g/V)*cos(myAC.mu_x + theta - alpha) ...
       -(g/V)*((rho*V^2)/(2*(myAC.W/myAC.S))) ...
          *( myAC.CL_alpha*alpha + myAC.CL_delta_e*delta_e(t) ...
          + myAC.CL_delta_s*delta_s(t)))...
             *1/( 1. + myAC.CL_alpha_dot*(myAC.mac/myAC.b)/(4*mu_rel));

% eq. (7.61c)
F3 = (myAC.Cm_0 + myAC.Cm_alpha*alpha ...
       + myAC.Cm_T_0 + myAC.Cm_T_alpha*alpha...      
       + myAC.Cm_delta_s*delta_s(t) ...
       + myAC.Cm_delta_e*delta_e(t) ...
       +(myAC.mac/(2*V))*myAC.Cm_q*q ...
       +(myAC.mac/(2*V))*myAC.Cm_alpha_dot*F2)...
          *((V/myAC.k_y)^2*(myAC.mac/myAC.b)/(2*mu_rel));

% eq. (7.61d)
F4 = V*cos(theta + myAC.mu_x - alpha);

% eq. (7.61e)
F5 = -V*sin(theta + myAC.mu_x - alpha);

% eq. (7.61f)
F6 = q;

% Costruzione della funzione di costo
f = [F1;
     F2;
     F3;
     F4;
     F5;
     F6];


end