%% FUNCTION USED TO FIND THE SADDLE POINT
% Use fsolve of the function to find the saddle point of the energy barrier
% theta_h = x(1), field angle for the hard layer
% theta_s = x(2), field angle for the soft layer
function out = direct_root(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor)
out = [sin(x(1))*cos(x(1)) + hp_h*sin(x(1)-theta_H_h) - j_exch*sin(x(2)-x(1)) + h_dip*(2*s2_factor*cos(x(2))*sin(x(1)) + s1_factor*sin(x(2)).*cos(x(1)));
    k*sin(x(2))*cos(x(2)) + m*hp_s*sin(x(2)-theta_H_s) + j_exch*sin(x(2)-x(1)) + h_dip*(2*s2_factor*sin(x(2))*cos(x(1))+ s1_factor*cos(x(2)).*sin(x(1)))];
end