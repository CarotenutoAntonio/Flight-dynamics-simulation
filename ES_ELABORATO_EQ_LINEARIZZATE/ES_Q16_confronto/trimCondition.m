function [alpha0_rad, alpha_0_deg, delta_e_0_rad, delta_e_0_deg, delta_s0_rad, delta_s_0_deg, delta_T0, theta0_rad, theta_0_deg, fval] = trimCondition(x0, lb, ub, gamma_0)
    options = optimset ('tolfun', 1e-9, 'Algorithm', 'interior-point');
    
    [x, fval] = fmincon(@costLongEquilibriumStaticStickFixed, ... % objectivefunc
        x0, ... initialguess
        [], [], [], [], ... linear constraints
        lb, ub, ... lower upper bounds
        @myNonLinearConstraint, ... nonlinear const.
        options); % mim serarch config options
    
    alpha0_rad = x(1);
    alpha_0_deg = convang(x(1), 'rad', 'deg');
    delta_e_0_rad = x(2);
    delta_e_0_deg = convang(x(2), 'rad', 'deg');
    delta_s0_rad = x(3);
    delta_s_0_deg = convang(x(3), 'rad', 'deg');
    delta_T0 = x(4);
    theta0_rad = gamma_0 - alpha0_rad;
    theta_0_deg = convang(theta0_rad, 'rad', 'deg');
end
