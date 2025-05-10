
function M = MassStickFree(t,x)
% termini che dipendono dalle variabili di stato
global g... 
       myAC 
   V           = x(1);
   alpha       = x(2);
   q           = x(3);
   x_EG        = x(4);
   z_EG        = x(5);
   theta       = x(6);
   delta_e_dot = x(7);
   delta_e     = x(8);
   
   rho = density(-z_EG);
   mu_rel = (myAC.W/g)/(rho*myAC.S*myAC.b);
   M = eye(8);
   
M(3,2) = -(1/(4*mu_rel)*(myAC.mac^2)/(myAC.k_y^2)*(V)/(myAC.b))*myAC.Cm_alpha_dot;
   
M(3,8) = -(1/(4*mu_rel))*(myAC.mac^2/myAC.k_y^2)*(V/myAC.b)*myAC.Cm_delta_e_dot;

M(7,1) = (myAC.mac_e*myAC.ec_adim)/myAC.k_e^2*sin(alpha-myAC.mu_x);

M(7,2) = ((myAC.mac_e*myAC.ec_adim)/myAC.k_e^2)*V*cos(alpha-myAC.mu_x)-...
         rho*V*myAC.S_e*myAC.mac_e^2/(4*myAC.I_e)*(1 - myAC.DepsDalpha)*myAC.Ch_e_alpha_dot;

M(7,3) = -(myAC.mac_e*myAC.ec_adim*myAC.x_C_e/myAC.k_e^2 - cos(myAC.Lambda_e));

M(7,6) = -(myAC.mac_e*myAC.ec_adim)/myAC.k_e^2*V*cos(alpha-myAC.mu_x);

M(7,8) = -rho*V*myAC.S_e*myAC.mac_e^2/(4*myAC.I_e)*myAC.Ch_e_delta_e_dot;

end
   

