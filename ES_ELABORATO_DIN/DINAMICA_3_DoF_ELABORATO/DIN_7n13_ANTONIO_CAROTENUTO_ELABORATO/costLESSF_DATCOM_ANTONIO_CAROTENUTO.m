function [f] = costLESSF_DATCOM_ANTONIO_CAROTENUTO(x)

%  Dichiarazione delle variabili globali
global g ...                          
       zEG_0 V0 q0 gamma0 ...        
       rho0 ...                      
       myAC                          

global...
    Adata WBflap v_alpha v_delta_e v_delta_flap v_Mach... 
Alpha_M Mach_a Deltae_M Mach_de Deltaf_M Mach_df... 
Alpha_de_M Deltae_a_M Mach_a_de... 
Alpha_df_M Deltaf_a_M Mach_a_df 


%Proprieta % dell ’aria alla quota attuale 
[air_Temp, sound_speed, air_pressure, rho0] = atmosisa(-zEG_0); 

%Densita % relativa del velivolo 
mu_rel = (myAC.W/g)/(rho0*myAC.S*myAC.b); 

%Il numero di Mach 
M0 = V0/sound_speed; 





%  Assegnazione delle componenti del design vector
alpha = x(1); 
delta_e = x(2); 
delta_flap = x(3); 
delta_T = x(4);

%Calcolo dei coefficienti aerodinamici 

[CL, CM, CD] = calcola_CL_CM_CD_ANTONIO_CAROTENUTO(alpha, delta_e, delta_flap,M0);
[Clq, Cmq] = calcola_CLq_CMq_ANTONIO_CAROTENUTO(M0);

%Riscrittura delle equazioni del moto , per volo stazionario 


F1 = (delta_T*myAC.T/myAC.W)*cos(alpha - myAC.mu_x + myAC.mu_T) ...
     -sin(gamma0)...
     -((rho0*V0^2)/(2*(myAC.W/myAC.S)))*CD ;


F2 = q0*(1 - ((myAC.mac/myAC.b)/(4*mu_rel))*Clq) ...
     -(delta_T*myAC.T/myAC.W)*(g/V0)...
     *sin(alpha - myAC.mu_x + myAC.mu_T)+...
     (g/V0)*cos(gamma0)-...
     ((rho0*V0^2)/(2*(myAC.W/myAC.S)))*(g/V0)...
     *CL;


F3 = (myAC.Cm_T_0 + myAC.Cm_T_alpha*alpha)*delta_T+...
     CM+(myAC.mac/(2*V0))*Cmq*q0;
    
%  Funzione obiettivo
f = F1*F1 + F2*F2 + F3*F3;

end
