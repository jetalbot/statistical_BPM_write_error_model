function h1 = h1Intgrandcone(v, u, phi, phih, gapsize, polesize, rh, zh, a, beta, zo, t)
d = gapsize/2 - t - zh;
A1 = polesize/2 + rh;
A2 = polesize/2 - rh;
B1 = cos(phi).*cos(phih) + beta.*sin(phi).*sin(phih);
C1 = gapsize/2 - d - u;
C2 = gapsize/2 + d + u;
ro = a.*(1 - u./zo);
r = v.*ro;
h1 = v.*ro.^2.*(log(((A1 - B1.*r).^2 + C2.^2)./( (A2 + B1.*r).^2 + C2.^2 )) - log(((A1 - B1.*r).^2 + C1.^2)./( (A2 + B1.*r).^2 + C1.^2 )));
end

