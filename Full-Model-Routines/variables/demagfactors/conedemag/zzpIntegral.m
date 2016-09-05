function out = zzpIntegral(phik, x, al, bl, t, alpha)
beta = bl./al;
zo = t./(1-alpha);

Iphik = cos(phik).^2 + 1./beta.^2.*sin(phik).^2;
A=(al./zo).^2 + Iphik;
B=(al./zo).^2.*cos(x) + Iphik;
z = 1-t./zo;

% I1t1 = 6.*log(A + sqrt(2*(A -B) ).*sqrt(A) -B) -6*log(A + sqrt(A).*sqrt(A.*(1 - t./zo).^2 -2*B.*(1 - t./zo) + A) - B.*(1 - t./zo)).*(1 - t./zo).^3;
% I1t2 = -2 + 2*(1 - t./zo).^3 - 3*(A.^2 -3.*B.^2)./A.^2.*log(-B + A + sqrt(A).*sqrt(2*(A - B))) + 3*(A.^2 - 3*B.^2)./A.^2.*log(-B + A.*(1 -t./zo) + sqrt(A).*sqrt(A.*(1 - t./zo).^2 -2*B.*(1 - t./zo) + A) );
% I1t3 = 3*(3*B + A).*sqrt(2*(A-B))./A.^(3/2) - 3*(3*B + A.*(1 -t./zo)).*sqrt(A.*(1 - t./zo).^2 -2*B.*(1 -t./zo) + A)./A.^(3/2);
% I1=zo.*B/(18*A.^(3/2)).*( I1t1 +I1t2 + I1t3 )

% I2t1 = -3*( A.^2 -3*B.^2)./A.^2.*log(-B.*(1 -t./zo) + A + sqrt(A).*sqrt(A.*((1 - t./zo).^2 +1) -2*B.*(1 - t./zo))).*(1 -t./zo).^3
% I2t2 = 3*(A.^2 - 3*B.^2)./A.^2.*log( -B.*(1 - t./zo) + A.*(1 -t./zo) + sqrt(A).*sqrt(2*A.*(1 -t./zo).^2 - 2*B.*(1 -t./zo).^2 )).*(1 -t./zo).^3
% I2t3 = 3*( 3*B.*(1 -t./zo) + A ).*sqrt(A.*((1-t./zo).^2 +1) -2*B.*(1 -t./zo)).*(1 -t./zo)./A.^(3/2) -2
% I2t4 = -3*(3*B.*(1 -t./zo) +A.*(1 -t./zo) ).*sqrt(2*(A - B).*(1 -t./zo).^2 ).*(1 - t./zo)./A.^(3/2) + 2*( 1 - t./zo).^3
% I2t5 = 6*log( A.*(1 -t./zo) -B + sqrt(A).*sqrt(A.*((1 -t./zo).^2 +1) -2*B.*(1 -t./zo)))
% I2t6 = -6*(1 -t./zo).^3.*log( (A -B).*(1 -t./zo) + sqrt(A).*sqrt(2*(A-B).*(1 -t./zo).^2) )
% I2=-zo.*B/(18*A.^(3/2)).*( I2t1 +I2t2 + I2t3 +I2t4 +I2t5 + I2t6 )
% 
% I3t1 = sqrt(A).*sqrt(2*(A-B)).*(4*A.^2 -B.*A -3*B.^2) + 3*B.*(A.^2 - B.^2).*log(-B + A + sqrt(A).*sqrt(2*(A-B)))
% I3t2 = -sqrt(A).*sqrt(A.*(1 -t./zo).^2 -2*B.*(1 -t./zo) + A).*( 2*((1 - t./zo).^2 +1).*A.^2 -B.*A.*(1 -t./zo) -3*B.^2 )
% I3t3 = -3*B.*( A.^2 -B.^2).*log( -B + A.*(1 - t./zo) + sqrt(A).*sqrt(A.*(1 -t./zo).^2) -2*B.*(1 -t./zo) + A )
% I3 =  zo./(6*A.^(7/2)).*(I3t1 + I3t2 + I3t3)
% 
% I4t1 = 3*B.*(A.^2 -B.^2).*log(-B.*(1 -t./zo) + A + sqrt(A).*sqrt( A.*((1 -t./zo).^2 +1) - 2*B.*(1 -t./zo) )).*(1 -t./zo).^3
% I4t2 = -3*B.*(A.^2 -B.^2).*log((A-B).*(1 -t./zo) + sqrt(A).*sqrt( 2*(A -B).*(1 -t./zo).^2 )).*(1 -t./zo).^3
% I4t3 = sqrt(A).*sqrt(A.*((1 -t./zo).^2 + 1) -2*B.*(1 -t./zo) ).*( 2*((1 - t./zo).^2 +1).*A.^2 -B.*A.*(1 -t./zo) -3*B.^2.*(1 -t./zo).^2 )
% I4t4 = -sqrt(A).*sqrt( 2*(A -B).*(1 -t./zo).^2 ).*( 4*(1 - t./zo).^2.*A.^2 -B.*A.*(1 -t./zo).^2 -3*B.^2.*(1 -t./zo).^2 )
% I4 = -zo./(6*A.^(7/2)).*(I4t1 + I4t2 + I4t3 + I4t4)

I1t1 = 6.*log(A + sqrt(2*(A -B) ).*sqrt(A) -B) -6*log(A + sqrt(A).*sqrt(A.*z.^2 -2*B.*z + A) - B.*z).*z.^3;
I1t2 = -2 + 2*z.^3 - 3*(A.^2 -3.*B.^2)./A.^2.*log(-B + A + sqrt(A).*sqrt(2*(A - B))) + 3*(A.^2 - 3*B.^2)./A.^2.*log(-B + A.*z + sqrt(A).*sqrt(A.*z.^2 -2*B.*z + A) );
I1t3 = 3*(3*B + A).*sqrt(2*(A-B))./A.^(3/2) - 3*(3*B + A.*z).*sqrt(A.*z.^2 -2*B.*z + A)./A.^(3/2);
I1 = B./(18*A.^(3/2)).*( I1t1 +I1t2 + I1t3 );

I2t1 = -3*( A.^2 -3*B.^2)./A.^2.*log(-B.*z + A + sqrt(A).*sqrt(A.*(z.^2 +1) -2*B.*z)).*z.^3;
I2t2 = 3*(A.^2 - 3*B.^2)./A.^2.*log( -B.*z + A.*z + sqrt(A).*sqrt(2*A.*z.^2 - 2*B.*z.^2 )).*z.^3;
I2t3 = 3*( 3*B.*z + A ).*sqrt(A.*(z.^2 +1) -2*B.*z).*z./A.^(3/2) -2;
I2t4 = -3*(3*B.*z +A.*z ).*sqrt(2*(A - B).*z.^2 ).*z./A.^(3/2) + 2*z.^3;
I2t5 = 6*log( A.*z -B + sqrt(A).*sqrt(A.*(z.^2 +1) -2*B.*z));
I2t6 = -6*z.^3.*log( (A -B).*z + sqrt(A).*sqrt(2*(A-B).*z.^2) );
I2 = -B./(18*A.^(3/2)).*( I2t1 +I2t2 + I2t3 +I2t4 +I2t5 + I2t6 );

I3t1 = sqrt(A).*sqrt(2*(A-B)).*(4*A.^2 -B.*A -3*B.^2) + 3*B.*(A.^2 - B.^2).*log(-B + A + sqrt(A).*sqrt(2*(A-B)));
I3t2 = -sqrt(A).*sqrt(A.*z.^2 -2*B.*z + A).*( 2*(z.^2 +1).*A.^2 -B.*A.*z -3*B.^2 );
I3t3 = -3*B.*( A.^2 -B.^2).*log( -B + A.*z + sqrt(A).*sqrt(A.*z.^2 -2*B.*z + A) );
I3 =  1./(6*A.^(7/2)).*(I3t1 + I3t2 + I3t3);

I4t1 = 3*B.*(A.^2 -B.^2).*log(-B.*z + A + sqrt(A).*sqrt( A.*(z.^2 +1) - 2*B.*z )).*z.^3;
I4t2 = -3*B.*(A.^2 -B.^2).*log((A-B).*z + sqrt(A).*sqrt( 2*(A -B).*z.^2 )).*z.^3;
I4t3 = sqrt(A).*sqrt(A.*(z.^2 + 1) -2*B.*z ).*( 2*(z.^2 +1).*A.^2 -B.*A.*z -3*B.^2.*z.^2 );
I4t4 = -sqrt(A).*sqrt( 2*(A -B).*z.^2 ).*( 4*z.^2.*A.^2 -B.*A.*z.^2 -3*B.^2.*z.^2 );
I4 = -1./(6*A.^(7/2)).*(I4t1 + I4t2 + I4t3 + I4t4);

out = I1 + I2 + I3 + I4; 
end