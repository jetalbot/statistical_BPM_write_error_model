% Calculates the S_factors used for calculating magnetostatic interaction
% between two ECC layers, both factors depend on relative position and
% shape. A circular cylinder is assumed here. These factors are used in the
% energy equation as:

% $$E = K_h V_h\sin^2(\theta_h) + K_s V_s\sin^2(\theta_s) -
% JA\cos(\theta_s - \theta_h) - \mu_0 M_h V_h\cos(\theta_h - \theta_H)$$
%
% $$-\mu_0 M_s V_s\cos(\theta_s - \theta_H)-
% \mu_0\frac{M_s V_s M_h V_h}{4\pi d^3}[2 S_2\cos(\theta_s)cos(\theta_h) -
% S_1\cos(\theta_s)cos(\theta_h)]$$

% where S_1 and S_2 are factors in the last term. See Bellegia et al. [1]
% for more details on the derivation of the calculation. 

% Z = distance between moments in nm (including spacing)
% R_rho = relative position
% a_1 = semi-major axis of first layer (we assume circular rather than
% elliptical cylinder)
% a_2 = semi-major axis of second layer
% t_1 = thickness of first layer
% t_2 = thickness of second layer

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

function out = S_1_integrand(k_rho, a_1, a_2, t_1, t_2, R_rho, Z)
z_val = k_rho.*R_rho;
z_val_1 = k_rho.*a_1;
z_val_2 = k_rho.*a_2;

[J_1] = besselj(1,z_val);
[J_1_1] = besselj(1,z_val_1);
[J_1_2] = besselj(1,z_val_2);

% ierr_s1
% ierr1_s1
% ierr2_s1
k_z_intgral_num_1 = S_1_k_z_intgral_num(abs(Z + (t_1 + t_2)/2), k_rho);
k_z_intgral_num_2 = S_1_k_z_intgral_num(abs(Z + (t_1 - t_2)/2), k_rho);
k_z_intgral_num_3 = S_1_k_z_intgral_num(abs(Z - (t_1 - t_2)/2), k_rho);
k_z_intgral_num_4 = S_1_k_z_intgral_num(abs(Z - (t_1 + t_2)/2), k_rho);

k_z_intgral_num = -k_z_intgral_num_1 + k_z_intgral_num_2 + k_z_intgral_num_3 - k_z_intgral_num_4;

out = J_1.*J_1_1.*J_1_2.*k_z_intgral_num./(k_rho.^3);

end

function out = S_1_k_z_intgral_num(alpha_p, beta_p)
% out = pi*(1 - alpha_p.*beta_p - exp(-alpha_p.*beta_p))/8;
out = pi*(-alpha_p.*beta_p - exp(-alpha_p.*beta_p))/8;
end

function out = S_2_integrand(k_rho, a_1, a_2, t_1, t_2, R_rho, Z)

z_val = k_rho.*R_rho;
z_val_1 = k_rho.*a_1;
z_val_2 = k_rho.*a_2;

[J_0] = besselj(0,z_val);
[J_1_1] = besselj(1,z_val_1);
[J_1_2] = besselj(1,z_val_2);

% ierr_s2
% ierr1_s2
% ierr2_s2

% J_0
k_z_intgral_num_1 = S_2_k_z_intgral_num(abs(Z + (t_1 + t_2)/2), k_rho);
k_z_intgral_num_2 = S_2_k_z_intgral_num(abs(Z + (t_1 - t_2)/2), k_rho);
k_z_intgral_num_3 = S_2_k_z_intgral_num(abs(Z - (t_1 - t_2)/2), k_rho);
k_z_intgral_num_4 = S_2_k_z_intgral_num(abs(Z - (t_1 + t_2)/2), k_rho);

k_z_intgral_num = -k_z_intgral_num_1 + k_z_intgral_num_2 + k_z_intgral_num_3 - k_z_intgral_num_4;

out = J_0.*J_1_1.*J_1_2.*k_z_intgral_num./(k_rho.^2);

end

function out = S_2_k_z_intgral_num(alpha_p, beta_p)
out = pi*exp(-alpha_p.*beta_p)/8;
end

function out = S_3_integrand(k_rho, a_1, a_2, t_1, t_2, R_rho, Z)

z_val = k_rho.*R_rho;
z_val_1 = k_rho.*a_1;
z_val_2 = k_rho.*a_2;

[J_2] = besselj(2,z_val);
[J_1_1] = besselj(1,z_val_1);
[J_1_2] = besselj(1,z_val_2);

k_z_intgral_num_1 = S_1_k_z_intgral_num(abs(Z + (t_1 + t_2)/2), k_rho);
k_z_intgral_num_2 = S_1_k_z_intgral_num(abs(Z + (t_1 - t_2)/2), k_rho);
k_z_intgral_num_3 = S_1_k_z_intgral_num(abs(Z - (t_1 - t_2)/2), k_rho);
k_z_intgral_num_4 = S_1_k_z_intgral_num(abs(Z - (t_1 + t_2)/2), k_rho);

k_z_intgral_num = -k_z_intgral_num_1 + k_z_intgral_num_2 + k_z_intgral_num_3 - k_z_intgral_num_4;

out = J_2.*J_1_1.*J_1_2.*k_z_intgral_num./(k_rho.^2);

end

function out = S_4_integrand(k_rho, a_1, a_2, t_1, t_2, R_rho, Z)

z_val = k_rho.*R_rho;
z_val_1 = k_rho.*a_1;
z_val_2 = k_rho.*a_2;

[J_1] = besselj(1,z_val);
[J_1_1] = besselj(1,z_val_1);
[J_1_2] = besselj(1,z_val_2);

k_z_intgral_num_1 = sign(Z + (t_1 + t_2)/2).*S_4_k_z_intgral_num(abs(Z + (t_1 + t_2)/2), k_rho);
k_z_intgral_num_2 = sign(Z + (t_1 - t_2)/2).*S_4_k_z_intgral_num(abs(Z + (t_1 - t_2)/2), k_rho);
k_z_intgral_num_3 = sign(Z - (t_1 - t_2)/2).*S_4_k_z_intgral_num(abs(Z - (t_1 - t_2)/2), k_rho);
k_z_intgral_num_4 = sign(Z - (t_1 + t_2)/2).*S_4_k_z_intgral_num(abs(Z - (t_1 + t_2)/2), k_rho);

k_z_intgral_num = -k_z_intgral_num_1 + k_z_intgral_num_2 + k_z_intgral_num_3 - k_z_intgral_num_4;

out = J_1.*J_1_1.*J_1_2.*k_z_intgral_num./(k_rho.^2);

end

function out = S_4_k_z_intgral_num(alpha_p, beta_p)
out = pi*(1 - exp(-alpha_p.*beta_p))/8;
end

%%
% References 
% [1] Beleggia, M., S. Tandon, Y. Zhu, and M. De Graef. "On the
% magnetostatic interactions between nanoparticles of arbitrary shape."
% Journal of magnetism and magnetic materials 278, no. 1 (2004): 270-284.