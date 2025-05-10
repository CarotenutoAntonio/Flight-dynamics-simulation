function [C_pitch] = Cpitch(alpha_B,delta_e,q)

C_pitch = -0.1885 - 0.00437*alpha_B - 0.123*q  - 0.0196*delta_e;

end