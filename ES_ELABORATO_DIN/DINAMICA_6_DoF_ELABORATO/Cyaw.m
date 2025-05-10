function [C_yaw] = Cyaw(alpha_B,beta_der,delta_a,delta_r,r)

if (alpha_B >= -5) && (alpha_B <= 10)
   C_yaw0 = 0.00125*beta_der;
elseif (alpha_B > 10) && (alpha_B <= 25)
       C_yaw0 = (0.00342 - 0.00022*alpha_B)*beta_der;
elseif (alpha_B > 25) && (alpha_B <= 35)
       C_yaw0 = -0.00201*beta_der;
else
    disp(['L''angolo di attacco eccede i limiti. Terminazione.']);
    return
end

C_yaw = C_yaw0 - 0.0142*r + delta_a/25*(0.00128 + 0.00213*alpha_B) + delta_r/30*(-0.0474 + 0.000804*alpha_B);

end