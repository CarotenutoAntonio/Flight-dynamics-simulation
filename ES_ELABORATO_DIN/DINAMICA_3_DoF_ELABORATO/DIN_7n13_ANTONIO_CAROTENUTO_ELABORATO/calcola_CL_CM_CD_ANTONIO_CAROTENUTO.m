function [CL, CM, CD] = calcola_CL_CM_CD_ANTONIO_CAROTENUTO(A ,DE ,DF ,M) 
%Richiamo le variabili globali 
global... 
Adata ... 
WBflap ... 
Alpha_M ...
Mach_a ... 
Deltae_M ...
Mach_de ... 
Deltaf_M ...
Mach_df... 
Alpha_de_M...
Deltae_a_M...
Mach_a_de...
Alpha_df_M... 
Deltaf_a_M... 
Mach_a_df 



%Si osservi che Mach_a e Alpha_M sono i valori in cui conosciamo i 
%valori delle funzioni mentre M e A sono i valori dove le interpoliamo, ossia
%dove estraiamo i valori di interesse 


CL_clear = interp2(Mach_a, Alpha_M, Adata.cl, M, A, 'spline'); 

CD_clear = interp2(Mach_a, Alpha_M, Adata.cd, M, A, 'spline'); 

CM_clear = interp2(Mach_a, Alpha_M, Adata.cm, M, A, 'spline'); 


%Incrementi dovuti alla deflessione dell %equilibratore 

dCL_de = interp2(Mach_de, Deltae_M, Adata.dcl_sym, M, DE, 'spline'); 

dCM_de = interp2(Mach_de, Deltae_M, Adata.dcm_sym, M, DE, 'spline'); 

dCDmin_de = interp2(Mach_de, Deltae_M, Adata.dcdmin_sym, M, DE, 'spline'); 

dCDi_de = interpn(Alpha_de_M, Deltae_a_M, Mach_a_de, Adata.dcdi_sym, A, DE, M, 'spline'); 


%Incrementi dovuti alla deflessione dei flap 
dCL_df = interp2(Mach_df, Deltaf_M, WBflap.dcl_sym, M, DF, 'spline');
dCM_df = interp2(Mach_df, Deltaf_M, WBflap.dcm_sym, M, DF, 'spline');
dCDmin_df = interp2(Mach_df, Deltaf_M, WBflap.dcdmin_sym, M, DF, 'spline');
dCDi_df = interpn(Alpha_df_M, Deltaf_a_M, Mach_a_df,WBflap.dcdi_sym, A, DF, M, 'spline');


%Calcolo dei coefficienti 
CL = CL_clear + dCL_de + dCL_df;

CD = CD_clear + dCDi_de + dCDmin_de + dCDi_df + dCDmin_df; 

CM = CM_clear + dCM_de + dCM_df;  

end