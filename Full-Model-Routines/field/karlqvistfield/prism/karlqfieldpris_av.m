function [havx havy havz] = karlqfieldpris_av(hg, phih, gapsize, polesize, d, rh, a, b, t, tol)

[h1 h2] = hIntegralpris(phih, gapsize, polesize, d, rh, a, b, t, tol);
vol = 4*a*b*t;
havx = -hg.*cos(phih).*h1./(2*pi*vol);
havy = -hg.*sin(phih).*h1./(2*pi*vol);
havz = hg.*h2./(pi*vol);
end