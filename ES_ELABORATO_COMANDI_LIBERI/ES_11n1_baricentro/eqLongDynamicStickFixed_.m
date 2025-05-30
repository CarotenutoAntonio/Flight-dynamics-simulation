function [dydt] = eqLongDynamicStickFixed_(t,state)
global g...                   %Accelerazione di gravità 
       zEG_0 V0 q0 gamma0...  %Condizioni iniziali
       rho0 ...               %Densità dell'aria all'altitudine h = (-zEG_0)
        myAC...                   %Oggetto 'Velivolo'
        delta_e...
        delta_s...
        delta_T
 
  % componenti del vettore di stato 
  V = state(1);
  alpha = state(2);
  q     = state(3);
  x_eg = state(4);
  zEG = state(5);
  theta = state(6);
 
    

      
  
  % 3DOF EQUATIONS 
  
  %% Eq. (7.63a) 
T_over_W = myAC.T/myAC.W;
 
alphaB_plus_muT = alpha - myAC.mu_x + myAC.mu_T;

gamma_t = theta - alpha + myAC.mu_x;

rhoVsquare_over_TwoWeightSurface = density(-zEG)*V^2/(2*(myAC.W/myAC.S));

CL_coeff = myAC.CL_alpha*alpha + myAC.CL_delta_e*delta_e(t)+...
           myAC.CL_delta_s*delta_s(t);

          V_dot = g*((delta_T(t)*T_over_W)*cos(alphaB_plus_muT)-...
           sin(gamma_t) - (rhoVsquare_over_TwoWeightSurface)*...
                          (myAC.CD_0 + myAC.K*(CL_coeff)^myAC.m));
 
%% Eq. (7.63b)
mu_rel = (myAC.W/g)/(density(-zEG)*myAC.S*myAC.b);    

macOverb_over_FourMuRel_dot_CLAlphadot = (myAC.mac/myAC.b)*myAC.CL_alpha_dot/(4*mu_rel);

macOverb_over_FourMuRel_dot_CLq = (myAC.mac/myAC.b)*myAC.CL_q/(4*mu_rel);

g_over_V = g/V;

           alpha_dot = 1/(1 + macOverb_over_FourMuRel_dot_CLAlphadot)*...
            (q*(1 - macOverb_over_FourMuRel_dot_CLq)-...
               delta_T(t)*(T_over_W)*(g_over_V)*sin(alphaB_plus_muT)+...
               (g_over_V)*cos(gamma_t)-...
               (rhoVsquare_over_TwoWeightSurface*g_over_V)*CL_coeff);
       
%%  Eq. (7.63c)
V_over_radiusgirationSquared = (V/myAC.k_y)^2;

macOverb_over_TwoMuRel = (myAC.mac/myAC.b)/(2*mu_rel);

CM_T = (myAC.Cm_T_0 + myAC.Cm_T_alpha*alpha)*delta_T(t);

CM_coeff = myAC.Cm_0 + myAC.Cm_alpha*alpha + myAC.Cm_delta_e*delta_e(t)+...
           myAC.Cm_delta_s*delta_s(t)+...
           (myAC.mac/(2*V))*(myAC.Cm_alpha_dot*alpha_dot + myAC.Cm_q*q); 
       
      q_dot = (V_over_radiusgirationSquared)*(macOverb_over_TwoMuRel)*...
                  (CM_T + CM_coeff);
       
% Eq. (7.63d)   
xEG_dot = V*cos(gamma_t);

% Eq. (7.63e)
zEG_dot = -V*sin(gamma_t);

% Eq. (7.63f)
theta_dot = q;
      
dydt = [V_dot;...
           alpha_dot;...
           q_dot;...
           xEG_dot;...
           zEG_dot;...
           theta_dot];
end 
function [rho] = density(h) 

[T,a,p,rho] = atmosisa(h); 

end 
