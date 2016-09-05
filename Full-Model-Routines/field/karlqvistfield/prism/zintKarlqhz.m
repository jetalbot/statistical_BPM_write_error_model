function intgral = zintKarlqhz(z, A,B1,B2,C)
intgral1a = z.*atan(A./(B1 - z)) + B1.*atan((B1 - z)./A) - (A.*log(A.^2 + (B1 - z).^2))/2;
intgral1b = z.*atan(C./(B1 - z)) + B1.*atan((B1 - z)./C) - (C.*log(C.^2 + (B1 - z).^2))/2;
intgral2a = z.*atan(A./(B2 + z)) - B2.*atan((B2 + z)./A) + (A.*log(A.^2 + (B2 + z).^2))/2;
intgral2b = z.*atan(C./(B2 + z)) - B2.*atan((B2 + z)./C) + (C.*log(C.^2 + (B2 + z).^2))/2;

intgral = intgral1a + intgral1b + intgral2a + intgral2b;
end