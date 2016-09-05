function intgral = karlqvistfieldrho(x, y, z, phih, g, G, rh, d)

A = g/2 + rh - (x.*cos(phih) + y.*sin(phih));
C = g/2 - rh + (x.*cos(phih) + y.*sin(phih));
B1 = G/2 - d;
B2 = G/2 + d;
% integral of log((A^2 + (B1 - x)^2)/(C^2 + (B1 - x)^2)) is 

% From mathematica: Integrate[Log[(a^2 + (b1 - z)^2)/(c^2 + (b1 - z)^2)], z] ==
% -2*a*ArcTan[(b1 - z)/a] + 2*c*ArcTan[(b1 - z)/c] - (b1 - z)*Log[(a^2 + (b1 - z)^2)/(c^2 + (b1 - z)^2)]

% integral of log((A^2 + (B2 + z)^2)/(C^2 + (B2 + z)^2)) is 
% From mathematica: Integrate[Log[(a^2 + (b2 + z)^2)/(c^2 + (b2 + z)^2)], z] ==
% 2*a*ArcTan[(b2 + z)/a] - 2*c*ArcTan[(b2 + z)/c] + (b2 + z)*Log[(a^2 + (b2 + z)^2)/(c^2 + (b2 + z)^2)]

% intgral1 = -2*A.*atan((B1 - z)./A) + 2*C.*atan((B1 - z)./C) - (B1 - z).*log((A.^2 + (B1 - z).^2)./(C.^2 + (B1 - z).^2));
% intgral2 =  2*A.*atan((B2 + z)./A) - 2*C.*atan((B2 + z)./C) + (B2 + z).*log((A.^2 + (B2 + z).^2)./(C.^2 + (B2 + z).^2));
% intgral = -intgral1 + intgral2;

intgral = zintKarlqhrho(z, A,B1,B2,C) - zintKarlqhrho(0, A,B1,B2,C);
end
