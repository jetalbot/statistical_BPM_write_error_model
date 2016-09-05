%% DETERMINANT OF THE ENERGY BARRIER HESSIAN MATRIX
% Returns the determinant of the energy barrier hessian matrix
function [delta e_x e_y e_xx e_yy e_xy] = hessian_det(x, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor)
% x(1) = th_s(1); % first minimum of the magnetisation angle theta
% x(2) = th_s(2); % second minimum of the magnetisation angle theta
% x(3) = hp_h; % applied field in the hard layer
% x(4) = hp_s; % applied field in the soft layer
% e_x = sin(x(1))*cos(x(1)) + x(3)*sin(x(1)-theta_H_h) - j_exch*sin(x(2)-x(1)) + h_dip*(2*s2_factor*cos(x(2))*sin(x(1)) + s1_factor*sin(x(2)).*cos(x(1)));
% e_y = k*sin(x(2))*cos(x(2)) + m*x(4)*sin(x(2)-theta_H_s)+ j_exch*sin(x(2)-x(1))+ h_dip*(2*s2_factor*sin(x(2))*cos(x(1)) + s1_factor*cos(x(2)).*sin(x(1)));


e_x = sin(x(1))*cos(x(1)) + x(3)*sin(x(1)-theta_H_h) - j_exch*sin(x(2)-x(1)) + h_dip*(2*s2_factor*cos(x(2))*sin(x(1)) + s1_factor*sin(x(2)).*cos(x(1)));
e_y = k*sin(x(2))*cos(x(2)) + m*x(4)*sin(x(2)-theta_H_s)+ j_exch*sin(x(2)-x(1))+ h_dip*(2*s2_factor*sin(x(2))*cos(x(1)) + s1_factor*cos(x(2)).*sin(x(1)));

e_xx = cos(x(1))^2 - sin(x(1))^2 + x(3)*cos(x(1) - theta_H_h) + j_exch*cos(x(2)-x(1)) + h_dip*(2*s2_factor*cos(x(2))*cos(x(1)) - s1_factor*sin(x(2)).*sin(x(1)));
e_yy = k*(cos(x(2))^2 - sin(x(2))^2) + m*x(4)*cos(x(2)- theta_H_s) + j_exch*cos(x(2)-x(1)) + h_dip*(2*s2_factor*cos(x(2))*cos(x(1))- s1_factor*sin(x(2)).*sin(x(1)));
e_xy  = -j_exch*cos(x(2) - x(1)) + h_dip*(-2*s2_factor*sin(x(2))*sin(x(1)) + s1_factor*cos(x(2))*cos(x(1)));
delta = e_xx*e_yy - e_xy*e_xy;
end