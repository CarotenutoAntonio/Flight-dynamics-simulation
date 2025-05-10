clc;
clear;
close all;

alldata=datcomimport('CITATION_SCRITTO_DA_ME_SENZA_SUP_MOBILI.out');

coeff_caso1=alldata{1,1};
coeff_caso2=alldata{1,2};
coeff_caso3=alldata{1,3};
coeff_caso4=alldata{1,4};

N=50;

%% CL
vcl_caso1=coeff_caso1.cl;
vcl_caso2=coeff_caso2.cl;
vcl_caso3=coeff_caso3.cl;
vcl_caso4=coeff_caso4.cl;

valpha=coeff_caso1.alpha;
valpha_graf=linspace(min(valpha),max(valpha),N);

vcl_caso1_graf=interp1(valpha,vcl_caso1,valpha_graf);
vcl_caso2_graf=interp1(valpha,vcl_caso2,valpha_graf);
vcl_caso3_graf=interp1(valpha,vcl_caso3,valpha_graf);
vcl_caso4_graf=interp1(valpha,vcl_caso4,valpha_graf);




figure(1)
plot(valpha_graf,vcl_caso1_graf,'k-',valpha_graf,vcl_caso2_graf,'y--',...
    valpha_graf,vcl_caso3_graf,'b-.',valpha_graf,vcl_caso4_graf,'r');

hold on
grid on
xlabel('$\alpha(deg)$','Interpreter','latex');
ylabel('$cl$','Interpreter','latex');
legend('wing','wing+body','wing+body+verticaltail','total');



%% CL_ALPHA
vcla_caso1=coeff_caso1.cla;
vcla_caso2=coeff_caso2.cla;
vcla_caso3=coeff_caso3.cla;
vcla_caso4=coeff_caso4.cla;

valpha=coeff_caso1.alpha;
valpha_graf=linspace(min(valpha),max(valpha),N);

vcla_caso1_graf=interp1(valpha,vcla_caso1,valpha_graf);
vcla_caso2_graf=interp1(valpha,vcla_caso2,valpha_graf);
vcla_caso3_graf=interp1(valpha,vcla_caso3,valpha_graf);
vcla_caso4_graf=interp1(valpha,vcla_caso4,valpha_graf);

figure(2)
plot(valpha_graf,vcla_caso1_graf,'k-',valpha_graf,vcla_caso2_graf,'y--',...
    valpha_graf,vcla_caso3_graf,'b-.',valpha_graf,vcla_caso4_graf,'r');

hold on
grid on
xlabel('$\alpha(deg)$','Interpreter','latex');
ylabel('$C_{L,\alpha}$','Interpreter','latex');
legend('wing','wing+body','wing+body+verticaltail','total');

%% CD

vcd_caso1=coeff_caso1.cd;
vcd_caso2=coeff_caso2.cd;
vcd_caso3=coeff_caso3.cd;
vcd_caso4=coeff_caso4.cd;

valpha=coeff_caso1.alpha;
valpha_graf=linspace(min(valpha),max(valpha),N);

vcd_caso1_graf=interp1(valpha,vcd_caso1,valpha_graf);
vcd_caso2_graf=interp1(valpha,vcd_caso2,valpha_graf);
vcd_caso3_graf=interp1(valpha,vcd_caso3,valpha_graf);
vcd_caso4_graf=interp1(valpha,vcd_caso4,valpha_graf);




figure(3)
plot(valpha_graf,vcd_caso1_graf,'k-',valpha_graf,vcd_caso2_graf,'g--',...
    valpha_graf,vcd_caso3_graf,'b-.',valpha_graf,vcd_caso4_graf,'r');

hold on
grid on
xlabel('$\alpha(deg)$','Interpreter','latex');
ylabel('$cd$','Interpreter','latex');
legend('wing','wing+body','wing+body+verticaltail','total');

%% CM
vcm_caso1=coeff_caso1.cm;
vcm_caso2=coeff_caso2.cm;
vcm_caso3=coeff_caso3.cm;
vcm_caso4=coeff_caso4.cm;

valpha=coeff_caso1.alpha;
valpha_graf=linspace(min(valpha),max(valpha),N);

vcm_caso1_graf=interp1(valpha,vcm_caso1,valpha_graf);
vcm_caso2_graf=interp1(valpha,vcm_caso2,valpha_graf);
vcm_caso3_graf=interp1(valpha,vcm_caso3,valpha_graf);
vcm_caso4_graf=interp1(valpha,vcm_caso4,valpha_graf);

figure(4)
plot(valpha_graf,vcm_caso1_graf,'k-',valpha_graf,vcm_caso2_graf,'y--',...
    valpha_graf,vcm_caso3_graf,'b-.',valpha_graf,vcm_caso4_graf,'r');

hold on
grid on
xlabel('$\alpha(deg)$','Interpreter','latex');
ylabel('$C_{M}$','Interpreter','latex');
legend('wing','wing+body','wing+body+verticaltail','total');


%% CM_ALPHA
vcma_caso1=coeff_caso1.cma;
vcma_caso2=coeff_caso2.cma;
vcma_caso3=coeff_caso3.cma;
vcma_caso4=coeff_caso4.cma;

valpha=coeff_caso1.alpha;
valpha_graf=linspace(min(valpha),max(valpha),N);

vcma_caso1_graf=interp1(valpha,vcma_caso1,valpha_graf);
vcma_caso2_graf=interp1(valpha,vcma_caso2,valpha_graf);
vcma_caso3_graf=interp1(valpha,vcma_caso3,valpha_graf);
vcma_caso4_graf=interp1(valpha,vcma_caso4,valpha_graf);

figure(5)
plot(valpha_graf,vcma_caso1_graf,'k-',valpha_graf,vcma_caso2_graf,'y--',...
    valpha_graf,vcma_caso3_graf,'b-.',valpha_graf,vcma_caso4_graf,'r');

hold on
grid on
xlabel('$\alpha(deg)$','Interpreter','latex');
ylabel('$C_{M,\alpha}$','Interpreter','latex');
legend('wing','wing+body','wing+body+verticaltail','total');

