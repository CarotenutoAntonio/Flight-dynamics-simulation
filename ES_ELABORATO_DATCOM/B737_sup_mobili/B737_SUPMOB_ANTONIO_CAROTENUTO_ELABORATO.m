% Importo i dati da Datcom mediante la funzione datcomimport 
 alldata1 = datcomimport('B737_ANTONIO_CAROTENUTO_ES_9_1_FLAP_rifatto.out', true, 0); 
 alldata2 = datcomimport('B737_ANTONIO_CAROTENUTO_ES_9_1_ELEVATOR_rifatto.out', true, 0); 

 caso1 = alldata1{1}; %flap 
 caso2 = alldata1{2}; %velivolo completo 
 
 caso11 = alldata2{1}; %flap 
 caso12 = alldata2{2}; %velivolo completo
%%%%%%%%%%%%%%%%%%%%%% CASO ELEVATORE FISSO %%%%%%%%%%%%%%%%%%%%%%%%% 

% Effetto sulla curva di portanza 
Df0_cl = caso1.dcl_sym(1); %effetto deflessione flap 0° 
Df20_cl = caso1.dcl_sym(2); %effetto deflessione flap 20° 
Df40_cl = caso1.dcl_sym(3); %effetto deflessione flap 40° 
De10_cl = caso2.dcl_sym; %effetto deflessione elevatore 10° 

% Effetto sulla curva di resistenza 
Df0_cd = caso1.dcdmin_sym(1); %effetto deflessione flap 0° 
Df20_cd = caso1.dcdmin_sym(2); %effetto deflessione flap 20° 
Df40_cd = caso1.dcdmin_sym(3); %effetto deflessione flap 40° 
De10_cd = caso2.dcdmin_sym; %effetto deflessione elevatore 10° 

% Effetto sulla curva di momento 
Df0_cm = caso1.dcm_sym(1); %effetto deflessione flap 0° 
Df20_cm = caso1.dcm_sym(2); %effetto deflessione flap 20° 
Df40_cm = caso1.dcm_sym(3); %effetto deflessione flap 40° 
De10_cm = caso2.dcm_sym; %effetto deflessione elevatore 10° 

%Effetto dovuto alla resistenza indotta dai flap e dall'equilibratore 
Df0_cd_indotta = caso1.dcdi_sym(:,1); %effetto deflessione flap 0°
Df20_cd_indotta = caso1.dcdi_sym(:,2); %effetto deflessione flap 20° 
Df40_cd_indotta = caso1.dcdi_sym(:,3); %effetto deflessione flap 40° 
De10_cd_indotta = caso2.dcdi_sym; %effetto deflessione elevatore 10° 

%I valori appena ottenuti, sommati alle curva del cl, cd, cm 
%restituiscono la curva di portanza, di resistenza e di momento in 
%presenza di tali deflessioni. 
% N.B. quanto appena detto vale nel campo lineare

%%%%%%%%%%%%%%%%%%%%%%%% PLOT CL-APLHA %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure; 
plot (caso2.alpha, caso2.cl+Df0_cl+De10_cl,'-b',... 
    caso2.alpha, caso2.cl+Df20_cl+De10_cl,'--r',... 
    caso2.alpha, caso2.cl+Df40_cl+De10_cl,'-g');

grid on; 
ylabel('C_L'); 
xlabel('\alpha (deg)'); 
title(['Curva di portanza (Mach = ' num2str(caso1.mach(1)) ')']); 
legend({'\delta_f = 0°', '\delta_f = 20°', '\delta_f = 40°'},... 
    'Location', 'southeast'); 

%%%%%%%%%%%%%%%%%%%%%%%% PLOT CD-APLHA %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure; 
plot (caso2.alpha, caso2.cd+Df0_cd+De10_cd+Df0_cd_indotta+De10_cd_indotta,'.-b',... 
caso2.alpha, caso2.cd+Df20_cd+De10_cd+Df20_cd_indotta+De10_cd_indotta,'--r',... 
caso2.alpha, caso2.cd+Df40_cd+De10_cd+Df40_cd_indotta+De10_cd_indotta,'-g'); 

grid on; 
ylabel('C_D'); 
xlabel('\alpha (deg)'); 
title(['Curva di portanza (Mach = ' num2str(caso1.mach(1)) ')']); 
legend({'\delta_f = 0°', '\delta_f = 20°', '\delta_f = 40°'},... 
'Location', 'southeast'); 


%%%%%%%%%%%%%%%%%%%%%%%% PLOT CM-APLHA %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure; 
plot (caso2.alpha, caso2.cm+Df0_cm+De10_cm,'-b',... 
    caso2.alpha, caso2.cm+Df20_cm+De10_cm,'--r',... 
    caso2.alpha, caso2.cm+Df40_cm+De10_cm,'-g');

grid on; 
ylabel('C_M'); 
xlabel('\alpha (deg)'); 
title(['Curva di momento (Mach = ' num2str(caso1.mach(1)) ')']); 
legend({'\delta_f = 0°', '\delta_f = 20°', '\delta_f = 40°'},... 
    'Location', 'southeast'); 


%% CASO FLAP FISSO e vari angoli di deflessione dell'equilibratore

% Effetto sulla curva di portanza 
Dem20_cl = caso11.dcl_sym(1); %effetto deflessione elevatore -20° 
De0_cl = caso11.dcl_sym(2); %effetto deflessione elevatore 0° 
Dep20_cl = caso11.dcl_sym(3); %effetto deflessione elevatore 20° 
Df15_cl = caso12.dcl_sym; %effetto deflessione flap 15° 

% Effetto sulla curva di resistenza 
Dem20_cd = caso11.dcdmin_sym(1); %effetto deflessione elevator -20° 
De0_cd = caso11.dcdmin_sym(2); %effetto deflessione elevator 0° 
Dep20_cd = caso11.dcdmin_sym(3); %effetto deflessione elevator 20° 
Df15_cd = caso12.dcdmin_sym; %effetto deflessione flap 15° 

% Effetto sulla curva di momento 
Dem20_cm = caso11.dcm_sym(1); %effetto deflessione elevator -20° 
De0_cm = caso11.dcm_sym(2); %effetto deflessione elevator 0° 
Dep20_cm = caso11.dcm_sym(3); %effetto deflessione elevator 20° 
Df15_cm = caso12.dcm_sym; %effetto deflessione flap 15° 

%Effetto dovuto alla resistenza indotta dai flap e dall'equilibratore 
Dem20_cd_indotta = caso11.dcdi_sym(:,1); %effetto deflessione elevator -20°
De0_cd_indotta = caso11.dcdi_sym(:,2); %effetto deflessione elevator 0° 
Dep20_cd_indotta = caso11.dcdi_sym(:,3); %effetto deflessione elevator 20° 
Df15_cd_indotta = caso12.dcdi_sym; %effetto deflessione elevator 15° 

%I valori appena ottenuti, sommati alle curva del cl, cd, cm 
%restituiscono la curva di portanza, di resistenza e di momento in 
%presenza di tali deflessioni. 
% N.B. quanto appena detto vale nel campo lineare

%%%%%%%%%%%%%%%%%%%%%%%% PLOT CL-APLHA %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure; 
plot (caso12.alpha, caso12.cl+Dem20_cl+Df15_cl,'-b',... 
    caso12.alpha, caso12.cl+De0_cl+Df15_cl,'--r',... 
    caso12.alpha, caso12.cl+Dep20_cl+Df15_cl,'-g');

grid on; 
ylabel('C_L'); 
xlabel('\alpha (deg)'); 
title(['Curva di portanza (Mach = ' num2str(caso11.mach(1)) ')']); 
legend({'\delta_e = -20°', '\delta_e = 0°', '\delta_e = -20°'},... 
    'Location', 'southeast'); 

%%%%%%%%%%%%%%%%%%%%%%%% PLOT CD-APLHA %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure; 
plot (caso12.alpha, caso12.cd+Dem20_cd+Df15_cd+Dem20_cd_indotta+Df15_cd_indotta,'.-b',... 
caso12.alpha, caso12.cd+De0_cd+Df15_cd+De0_cd_indotta+Df15_cd_indotta,'--r',... 
caso12.alpha, caso12.cd+Dep20_cd+Df15_cd+Dep20_cd_indotta+Df15_cd_indotta,'-g'); 

grid on; 
ylabel('C_D'); 
xlabel('\alpha (deg)'); 
title(['Curva di portanza (Mach = ' num2str(caso11.mach(1)) ')']); 
legend({'\delta_e = -20°', '\delta_e = 0°', '\delta_e = 20°'},... 
'Location', 'southeast'); 


%%%%%%%%%%%%%%%%%%%%%%%% PLOT CM-APLHA %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure; 
plot (caso12.alpha, caso12.cm+Dem20_cm+Df15_cm,'-b',... 
    caso12.alpha, caso12.cm+De0_cm+Df15_cm,'--r',... 
    caso12.alpha, caso12.cm+Dep20_cm+Df15_cm,'-g');

grid on; 
ylabel('C_M'); 
xlabel('\alpha (deg)'); 
title(['Curva di momento (Mach = ' num2str(caso11.mach(1)) ')']); 
legend({'\delta_e = -20°', '\delta_e = 0°', '\delta_e = 20°'},... 
    'Location', 'southeast'); 
