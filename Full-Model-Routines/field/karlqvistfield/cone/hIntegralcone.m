function [h1 h2] = hIntegralcone(phih, gapsize, polesize, rh, zh, a, beta, t, alpha, tol)
zo = t./(1-alpha);
phimin = 0; phimax = 2*pi;

umin = 0; umax = t;
vmin = 0; vmax = 1;

h1 = triplequad(@h1Intgrandcone, vmin, vmax, umin, umax, phimin, phimax, tol, [], phih, gapsize, polesize, rh, zh, a, beta, zo, t);
h2 = triplequad(@h2Intgrandcone, vmin, vmax, umin, umax, phimin, phimax, tol, [], phih, gapsize, polesize, rh, zh, a, beta, zo, t);
end