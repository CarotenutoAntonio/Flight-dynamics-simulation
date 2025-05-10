function [Croll] = Croll(alpha_B,beta_der,delta_a,delta_r,p,r)

if (alpha_B >= -5) && (alpha_B <= 15)
   Croll_0 = (-0.00092 - 0.00012*alpha_B)*beta_der;
elseif (alpha_B > 15) && (alpha_B <= 25)
       Croll_0 = (-0.006 + 0.00022*alpha_B)*beta_der;
else
    disp(['L''angolo di attacco eccede i limiti. Terminazione.']);
    return
end

Croll = Croll_0 - 0.315*p + 0.0126*r + (delta_a/25)*(-0.0628 + 0.00121*alpha_B) - (delta_r/30)*(-0.00124 + 0.000351*alpha_B);

end
