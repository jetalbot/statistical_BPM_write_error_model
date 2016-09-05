function out = S_1_k_z_intgral_num(alpha_p, beta_p)
% out = pi*(1 - alpha_p.*beta_p - exp(-alpha_p.*beta_p))/8;
out = pi*(-alpha_p.*beta_p - exp(-alpha_p.*beta_p))/8;
end