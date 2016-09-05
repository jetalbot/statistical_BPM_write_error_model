function [ nxx nyy nzz ] = conedemagIntegral(a, b, t, alpha, tol)
beta = b./a;
V= pi*t.*a.*b.*((alpha + 1/2).^2 + 3/4)/3;
zo = t./(1-alpha);
cxx = a.*b.*zo./(2*pi*V);
cyy = a.*b.*zo./(2*pi*V.*beta.^2);
czz = a.*b.*zo./(2*pi*V);

phikmin = 0; phikmax = 2*pi;
xmax = pi; xmin = 0.0000005; % 0

for j =1:length(a)
%     Ixx(j) = dblquad(@integrndxx,phikmin,phikmax,xmin,xmax,tolx,@myquad,a(j),b(j),t(j),alpha(j));
%     Iyy(j) = dblquad(@integrndyy,phikmin,phikmax,xmin,xmax,toly,@myquad,a(j),b(j),t(j),alpha(j));
%     Izz(j) = dblquad(@integrndzz,phikmin,phikmax,xmin,xmax,tolz,@myquad,a(j),b(j),t(j),alpha(j));

    Ixx(j) = dblquad(@integrndconexx,phikmin,phikmax,xmin,xmax,tol,[],a(j),b(j),t(j),alpha(j));
    Iyy(j) = dblquad(@integrndconeyy,phikmin,phikmax,xmin,xmax,tol,[],a(j),b(j),t(j),alpha(j));
    Izz(j) = dblquad(@integrndconezz,phikmin,phikmax,xmin,xmax,tol,[],a(j),b(j),t(j),alpha(j));

    nxxp(j) = cxx(j).*Ixx(j);
    nyyp(j) = cyy(j).*Iyy(j);
    nzzp(j) = czz(j).*Izz(j);
end
nzz = 1 -nzzp;
nyy = nyyp;
nxx = nxxp;