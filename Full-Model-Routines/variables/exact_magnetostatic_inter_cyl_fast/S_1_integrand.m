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