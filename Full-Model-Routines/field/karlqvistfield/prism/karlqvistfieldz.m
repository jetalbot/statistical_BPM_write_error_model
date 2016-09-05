function intgral = karlqvistfieldz(x, y, z, phih, g, G, rh, d)

A = g/2 + rh - (x.*cos(phih) + y.*sin(phih));
C = g/2 - rh + (x.*cos(phih) + y.*sin(phih));
B1 = G/2 - d;
B2 = G/2 + d;
% integral of atan(A/(B1 -z)) is 

% From mathematica: Integrate[ArcTan[a/(b1 - x)], x] ==
% x*ArcTan[a/(b1 - x)] + b1*ArcTan[(b1 - x)/a] - (a*Log[a^2 + (b1 - x)^2])/2

% Integrate[ArcTan[c/(b1 - x)], x] ==
% x*ArcTan[c/(b1 - x)] + b1*ArcTan[(b1 - x)/c] - (c*Log[b1^2 + c^2 - 2*b1*x + x^2])/2

% Integrate[ArcTan[a/(b2 + x)], x] ==
% x*ArcTan[a/(b2 + x)] - b2*ArcTan[(b2 + x)/a] + (a*Log[a^2 + (b2 + x)^2])/2

% Integrate[ArcTan[c/(b2 + x)], x] ==
% x*ArcTan[c/(b2 + x)] - b2*ArcTan[(b2 + x)/c] + (c*Log[b2^2 + c^2 + 2*b2*x + x^2])/2

intgral = zintKarlqhz(z, A,B1,B2,C) - zintKarlqhz(0, A,B1,B2,C);
end
