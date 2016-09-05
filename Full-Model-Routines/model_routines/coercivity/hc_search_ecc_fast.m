function [h tprob xo x exitflag] = hc_search_ecc_fast(theta_H_h, theta_H_s, k_h, v_h, j_exch, h_dip, k, m, kb, temp, xo, s1_factor, s2_factor)
ho = 0; %1e-6;

options = optimset('Display','iter','TolFun',1e-50, 'TolX',1e-15);
% [h, f, exitflag] = fsolve(@(x)switchrateIntegrand_function_ecc(x, theta_H_h, theta_H_s, k_h, v_h, j_exch, h_dip, k, m, kb, temp, xo, s1_factor, s2_factor), ho, options);
[h, f, exitflag] = fzero(@(x)switchrateIntegrand_function_ecc(x, theta_H_h, theta_H_s, k_h, v_h, j_exch, h_dip, k, m, kb, temp, xo, s1_factor, s2_factor), ho, options);

x =f;
h_h = h;
h_s = h;
tprob =switchrateIntegrand_ecc(h_h, h_s, theta_H_h, theta_H_s, k_h, v_h, j_exch, h_dip, k, m, kb, temp, s1_factor, s2_factor);
end
