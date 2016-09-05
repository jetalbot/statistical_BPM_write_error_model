function [h1 h2] = hIntegralcyl(phih, gapsize, polesize, rh, zh, a, beta, t, tol)
phimin = 0; phimax = 2*pi;

umin = 0; umax = t;
vmin = 0; vmax = 1;

h1 = triplequad(@h1Intgrandcyl, vmin, vmax, umin, umax, phimin, phimax, tol, [], phih, gapsize, polesize, rh, zh, a, beta, t);
h2 = triplequad(@h2Intgrandcyl, vmin, vmax, umin, umax, phimin, phimax, tol, [], phih, gapsize, polesize, rh, zh, a, beta, t);

end