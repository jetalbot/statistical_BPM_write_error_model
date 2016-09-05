function [havx havy havz] = karlqfieldcone_av(hg, phih, gapsize, polesize, rh, zh, a, b, t, alpha, tol)
beta = b./a;
[h1 h2] = hIntegralcone(phih, gapsize, polesize, rh, zh, a, beta, t, alpha, tol);
vol= pi*t.*a.*b.*((alpha + 1/2).^2 + 3/4)/3;
havx = -beta.*hg.*cos(phih).*h1./(2*pi*vol);
havy = -beta.*hg.*sin(phih).*h1./(2*pi*vol);
havz = beta.*hg.*h2./(pi*vol);

end