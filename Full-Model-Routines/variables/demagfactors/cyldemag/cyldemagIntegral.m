function [ nxx nyy nzz ] = cyldemagIntegral(a, b, t, tol)
beta = b./a;
V= pi*t.*a.*b;

cxx = a.*b./(2*pi*V);
cyy = a.*b./(2*pi*V.*beta.^2);
czz = a.*b./(2*pi*V);

phikmin = 0; phikmax = 2*pi;

for j =1:length(a)
    Ixx(j) = quad(@integrndcylxx,phikmin,phikmax,tol,[],a(j),b(j),t(j));
    Iyy(j) = quad(@integrndcylyy,phikmin,phikmax,tol,[],a(j),b(j),t(j));
    Izz(j) = quad(@integrndcylzz,phikmin,phikmax,tol,[],a(j),b(j),t(j));

    nxxp(j) = cxx(j).*Ixx(j);
    nyyp(j) = cyy(j).*Iyy(j);
    nzzp(j) = czz(j).*Izz(j);
end
nzz = 1 -nzzp;
nyy = nyyp;
nxx = nxxp;