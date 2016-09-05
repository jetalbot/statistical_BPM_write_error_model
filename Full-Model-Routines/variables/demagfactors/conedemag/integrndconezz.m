function out = integrndconezz(phik, x, a, b, t, alpha)

beta = b./a;

for j = 1:length(phik)
zzpI(j) = zzpIntegral(phik(j), x, a, b, t, alpha);
end
 
out = cos(x).*sqrt(cos(phik).^2 + sin(phik).^2./beta.^2 ).*zzpI;
end
