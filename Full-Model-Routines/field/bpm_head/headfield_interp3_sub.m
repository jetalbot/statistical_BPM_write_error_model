function h_interp = headfield_interp3_sub(x, y, z, nx, ny, nz, stepx, stepy, stepz, x_vect, y_vect, z_vect, h)
min_nx = max(1, min(nx - stepx(1)));
max_nx = min(length(x_vect), max(abs(nx) + stepx(2)));

min_ny = max(1, min(ny - stepy(1)));
max_ny = min(length(y_vect), max(abs(ny) + stepy(2)));

min_nz = max(1, min(nz - stepz(1)));
max_nz = min(length(z_vect), max(abs(nz) + stepz(2)));

xbit=x_vect(min_nx:max_nx);
ybit=y_vect(min_ny:max_ny);
zbit=z_vect(min_nz:max_nz);

hbit=h(min_ny:max_ny, min_nx:max_nx, min_nz:max_nz)

h_interp = interp3(xbit, ybit, zbit, hbit, x, y, z, 'cubic');
end
