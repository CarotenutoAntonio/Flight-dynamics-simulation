function [C_D] = CD(alpha_B)

if (alpha_B >= -5) && (alpha_B <= 20)
   C_D = 0.1423 - 0.00438*alpha_B + 0.0013*alpha_B^2;
elseif (alpha_B > 20) && (alpha_B <= 90)
       C_D = -0.358 - 0.0473*alpha_B + 0.0000348*alpha_B^2;
else
    disp(['L''angolo di attacco eccede i limiti. Terminazione']);
    return
end

end
