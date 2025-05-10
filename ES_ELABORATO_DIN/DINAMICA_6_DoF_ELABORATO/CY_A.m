function [CY_A] = CY_A(alpha_B,beta_der,delta_a,delta_r)

CY_A = -0.0186*beta_der + (delta_a/25)*(0.039 - 0.00227*alpha_B)+...
       (delta_r/30)*(0.141 - 0.00265*alpha_B);

end