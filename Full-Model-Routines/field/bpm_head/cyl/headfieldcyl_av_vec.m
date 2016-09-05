function [hx_av hy_av hz_av] = headfieldcyl_av_vec(a, b, t, r_vec_center, x_vect, y_vect, z_vect, step_vect, tol, hx, hy, hz)

hx_av = headfieldcyl_av(a, b, t, r_vec_center, x_vect, y_vect, z_vect, step_vect, tol, hx);
hy_av = headfieldcyl_av(a, b, t, r_vec_center, x_vect, y_vect, z_vect, step_vect, tol, hy);
hz_av = headfieldcyl_av(a, b, t, r_vec_center, x_vect, y_vect, z_vect, step_vect, tol, hz);
end