function [havx havy havz] = karlqfieldcyl_av(hg, phih, gapsize, polesize, rh, zh, a, b, t, tol)
beta = b./a;
[h1 h2] = hIntegralcyl(phih, gapsize, polesize, rh, zh, a, beta, t, tol);
vol= pi*a.*b.*t;
havx = -beta.*hg.*cos(phih).*h1./(2*pi*vol);
havy = -beta.*hg.*sin(phih).*h1./(2*pi*vol);
havz = beta.*hg.*h2./(pi*vol);

end