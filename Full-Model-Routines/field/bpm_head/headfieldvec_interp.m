function [hx_interp hy_interp hz_interp] = headfieldvec_interp(x, y, nx, ny, stepx, stepy, x_vect, y_vect, hx, hy, hz)

% min_nx = (min(nx - stepx(1))<=0) + min(nx - stepx(1)).*(min(nx - stepx(1))>0); % ensures that min_nx >=1
% max_nx = length(x_vect).*(max(nx + stepx(2))>length(x_vect)) + max(nx + stepx(2)).*(max(nx + stepx(2))<=length(x_vect)); % ensures that max_nx<=length(x_vect)
% 
% min_ny = (min(ny - stepy(1))<=0) + min(ny - stepy(1)).*(min(ny - stepy(1))>0); % ensures that min_ny >=1
% max_ny = length(y_vect).*(max(ny + stepy(2))>length(y_vect)) + max(ny + stepy(2)).*(max(ny + stepy(2))<=length(y_vect)); % ensures that max_ny<=length(y_vect)

min_nx = max(1, min(nx - stepx(1)));
max_nx = min(length(x_vect), max(abs(nx) + stepx(2)));

min_ny = max(1, min(ny - stepy(1)));
max_ny = min(length(y_vect), max(abs(ny) + stepy(2)));
xbit=x_vect(min_nx:max_nx);
ybit=y_vect(min_ny:max_ny);
hxbit=hx(min_ny:max_ny,min_nx:max_nx);
hybit=hy(min_ny:max_ny,min_nx:max_nx);
hzbit=hz(min_ny:max_ny,min_nx:max_nx);

hx_interp = interp2(xbit, ybit, hxbit, x, y,'*spline');
hy_interp = interp2(xbit, ybit, hybit, x, y,'*spline');
hz_interp = interp2(xbit, ybit, hzbit, x, y,'*spline');
end
