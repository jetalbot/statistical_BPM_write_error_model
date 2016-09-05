%% ENERGY BARRIER MINIMUM FUNCTION
function out = energy_min(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor)
theta_h = x(1); % magnetisation field angle in the hard layer
theta_s = x(2); % magnetisation field angle in the soft layer
out = (1/2)*(sin(theta_h).^2 + k.*sin(theta_s).^2) - hp_h.*cos(theta_h - theta_H_h) - m.*hp_s.*cos(theta_s - theta_H_s) - j_exch.*cos(theta_s -theta_h) - h_dip.*(2*s2_factor*cos(theta_s).*cos(theta_h) - s1_factor.*sin(theta_s).*sin(theta_h));
end