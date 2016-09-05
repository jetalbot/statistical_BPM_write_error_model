function h2 = h2Intgrandcone(v, u, phi, phih, gapsize, polesize, rh, zh, a, beta, zo, t)
d = gapsize/2 - t - zh;
A1 = polesize/2 + rh;
A2 = polesize/2 - rh;
B1 = cos(phi).*cos(phih) + beta.*sin(phi).*sin(phih);
C1 = gapsize/2 - d - u;
C2 = gapsize/2 + d + u;
ro = a.*(1 - u./zo);
r = v.*ro;
h2 = v.*ro.^2.*(atan((A1 - B1.*r)./C1) + atan((A2 + B1.*r)./C1) + atan((A1 - B1.*r)./C2) + atan((A2 + B1.*r)./C2));
end

