function out = integrndcylxx(phik, a, b,t)
beta = b./a;
A = cos(phik).^2 + sin(phik).^2./beta.^2;
[elipk,elipe]=ellipke(4*a.^2./(A.*t.^2 +4*a.^2));
I = -8.*a./(3.*sqrt(A)) + t.*(4/3 -A.*t.^2./(3*a.^2)).*sqrt((A.*t.^2 + 4*a.^2)./(A.*t.^2)).*elipe + t.*(4/3 + A.*t.^2./(3*a.^2)).*sqrt(A.*t.^2./(A.*t.^2 + 4*a.^2)).*elipk;
out = cos(phik).^2./A.*I;
end
