function [h1 h2] = hIntegralpris(phih, gapsize, polesize, d, rh, a, b, t, tol)

xmin = -a;
xmax = a;
ymin = -b;
ymax  = b;

h1 = dblquad(@karlqvistfieldrho, xmin, xmax, ymin, ymax, tol, [], t, phih, polesize, gapsize, rh, d);
h2 = dblquad(@karlqvistfieldz, xmin, xmax, ymin, ymax, tol, [], t, phih, polesize, gapsize, rh, d);

end