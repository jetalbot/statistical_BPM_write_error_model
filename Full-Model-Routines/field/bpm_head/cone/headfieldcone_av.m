function h_av = headfieldcone_av(a, b, t, alpha, r_vec_center, x_vect, y_vect, z_vect, step_vect, tol, h)
phimin = 0; phimax = 2*pi;
umin = 0; umax = t;
vmin = 0; vmax = 1;

h_intgral = triplequad(@hIntgrandcone, vmin, vmax, umin, umax, phimin, phimax, tol, [],  a, b, t, alpha, r_vec_center, x_vect, y_vect, z_vect, step_vect, h);
vol= pi*t.*a.*b.*((alpha + 1/2).^2 + 3/4)/3;
h_av = h_intgral./vol;
end