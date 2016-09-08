function tprob =switchrateIntegrand_ecc(h_h, h_s, theta_H_h, theta_H_s, k_h, v_h, j_exch, h_dip, k, m, kb, temp)
% function tprob = switchrateIntegrand_ecc(h, theta_h, k_h, v_h, j_exch, h_dip, k, m, kb, temp)
% theta_h = theta_h.*ones(size(h));
% % [eb1 eb2] = scaledbarrier_ecc(h, theta_h, k_h, v_h, j_exch, h_dip, k, m, kb, temp);
% h_h = h;
% h_s = h;
% theta_H_h = theta_h;
% theta_H_s = theta_h;

 theta_H_h = theta_H_h.*ones(size(h_h));
 theta_H_s = theta_H_s.*ones(size(h_s));
 
[eb1 eb2] = scaledbarrier_ecc(h_h, h_s, theta_H_h, theta_H_s, k_h, v_h, j_exch, h_dip, k, m, kb, temp);
tprob = exp(-eb1)+ exp(-eb2);
end