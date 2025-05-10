function [C_L] = CL(alpha_B,delta_e)

if (alpha_B >= -5) && (alpha_B <= 10)
   C_L = 0.732 + 0.0751*alpha_B + 0.0144*delta_e;
elseif (alpha_B > 10) && (alpha_B <= 40) %%
       C_L = 0.569 + 0.106*alpha_B - 0.00148*alpha_B^2 + 0.0144*delta_e;
else
    disp(['L''angolo di attacco eccede i limiti. Terminazione.']);
    return
end

end