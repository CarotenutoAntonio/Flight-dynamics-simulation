function [CLq, CMq] = calcola_CLq_CMq_ANTONIO_CAROTENUTO(Mach) 
global v_Mach Adata 
CLq = interp1(v_Mach, Adata.clq(1,:), Mach, 'spline'); 
CMq = interp1(v_Mach, Adata.cmq(1,:), Mach, 'spline'); 

%Consideriamo solo la prima riga dei valori Adata.clq e Adata.cmq poichè
%gli altri valori sono pari a 9999. Sono, cioè, valori non ricavabili da 
%DATCOM 
end 