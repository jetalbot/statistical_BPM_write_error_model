function h_av = headfieldpris_av(a, b, t, r_vec_center, x_vect, y_vect, z_vect, step_vect, tol, h)
xmin = -a; xmax = a;
ymin = -b; ymax = b;
zmin = 0; zmax = t;

h_intgral = triplequad(@hIntgrandpris, xmin, xmax, ymin, ymax, zmin, zmax, tol, [], r_vec_center, x_vect, y_vect, z_vect, step_vect, h);
vol = 4*a.*b.*t;
h_av = h_intgral./vol;
end