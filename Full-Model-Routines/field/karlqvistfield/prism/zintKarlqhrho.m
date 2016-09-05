function intgral = zintKarlqhrho(z, A,B1,B2,C)
intgral1 = -2*A.*atan((B1 - z)./A) + 2*C.*atan((B1 - z)./C) - (B1 - z).*log((A.^2 + (B1 - z).^2)./(C.^2 + (B1 - z).^2));
intgral2 =  2*A.*atan((B2 + z)./A) - 2*C.*atan((B2 + z)./C) + (B2 + z).*log((A.^2 + (B2 + z).^2)./(C.^2 + (B2 + z).^2));
intgral = -intgral1 + intgral2;
end