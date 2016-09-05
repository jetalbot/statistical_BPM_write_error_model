function [S_1_factor, S_2_factor, S_3_factor, S_4_factor] = S_factors_fast(Z, R_rho, a_1, a_2, t_1,t_2)

k_rho_min = 0;  
k_rho_max = Inf; 

rel_tol = 1e-8;
abs_tol = 1e-12;
max_cnt = 1000;

 
S_1_intgral = quadgk(@(x)S_1_integrand(x, a_1, a_2, t_1, t_2, R_rho, Z), k_rho_min, k_rho_max,'RelTol', rel_tol,'AbsTol', abs_tol,'MaxIntervalCount', max_cnt);

S_2_intgral = quadgk(@(x)S_2_integrand(x, a_1, a_2, t_1, t_2, R_rho, Z), k_rho_min, k_rho_max,'RelTol', rel_tol,'AbsTol', abs_tol, 'MaxIntervalCount', max_cnt);

S_3_intgral = quadgk(@(x)S_3_integrand(x, a_1, a_2, t_1, t_2, R_rho, Z), k_rho_min, k_rho_max,'RelTol', rel_tol,'AbsTol', abs_tol, 'MaxIntervalCount', max_cnt);

S_4_intgral = quadgk(@(x)S_4_integrand(x, a_1, a_2, t_1, t_2, R_rho, Z), k_rho_min, k_rho_max,'RelTol', rel_tol,'AbsTol', abs_tol, 'MaxIntervalCount', max_cnt);

V_1 = pi*a_1^2*t_1;
V_2 = pi*a_2^2*t_2;

R = sqrt(Z.^2 + R_rho.^2);
S_1_factor = 8*pi*R^3*4*a_1*a_2*S_1_intgral/(R_rho*V_1*V_2);
S_2_factor = 8*pi*R^5*4*a_1*a_2*S_2_intgral/((R_rho^2 - 2*Z^2)*V_1*V_2);
S_3_factor = 8*pi*R^5*4*a_1*a_2*S_3_intgral/(3*R_rho^2*V_1*V_2);
S_4_factor = 8*pi*R^5*4*a_1*a_2*S_4_intgral/(3*Z*R_rho*V_1*V_2);

end