function [h tprob xo x exitflag] = hc_search_fast(theta_h, h_keff, mu_o, m_s, vol, kb, temp, xo)
ho = 0; %1e-6;

options = optimset('Display','iter','TolFun',1e-50, 'TolX',1e-15);
% [h, f, exitflag] = fsolve(@(x)switchrateIntegrand_function_ecc(x, theta_H_h, theta_H_s, k_h, v_h, j_exch, h_dip, k, m, kb, temp, xo, s1_factor, s2_factor), ho, options);
[h, f, exitflag] = fzero(@(x)switchrateIntegrand_function(x, theta_h, h_keff, mu_o, m_s, vol, kb, temp, xo), ho, options);

x =f;
tprob =switchrateIntegrand(h, theta_h, h_keff, mu_o, m_s, vol, kb, temp);
end
