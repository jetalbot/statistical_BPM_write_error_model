function intgrand = hIntgrandcone(v, u, phi, a, b, t, alpha, r_vec_center, x_vect, y_vect, z_vect, step_vect, h)
zo = t/(1-alpha);
r = v.*a.*(1 - u./zo);

x = r.*cos(phi);
y = (b./a).*r.*sin(phi);
z = u;

x1 = r_vec_center(1) + x;
y1 = r_vec_center(2) + y;
z1 = r_vec_center(3) + z;
z1 = z1.*ones(size(x1));

% spacing dx, dy, dz is not always 1 so divide by dx and dy, possible get these
% from a vector parameter

if all(all(step_vect)) % check whether all elements of step_vect are non zero
    dx = x_vect(2) - x_vect(1);
    dy = y_vect(2) - y_vect(1);
    dz = z_vect(2) - z_vect(1);

    nx1 = floor((x1-x_vect(1))/dx);
    ny1 = floor((y1-y_vect(1))/dy);
    nz1 = floor((z1-z_vect(1))/dz);

    stepx = step_vect(1,:);
    stepy = step_vect(2,:);
    stepz = step_vect(3, :);

    h_interp = headfield_interp3_sub(x1, y1, z1, nx1, ny1, nz1, stepx, stepy, stepz, x_vect, y_vect, z_vect, h);
else

    h_interp = headfield_interp3(x1, y1, z1, x_vect, y_vect, z_vect, h);
end
intgrand = (b./a).*(a.^2).*v.*((1 - u./zo).^2).*h_interp;
end

