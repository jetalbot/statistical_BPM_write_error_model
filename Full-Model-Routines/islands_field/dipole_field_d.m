function [dxx dxy dxz dyx dyy dyz dzx dzy dzz] = dipole_field_d(x,y,z)

% All values in SI units
% (x,y,z) is the location of the field point relative to the dipole in m

% (x,y,z) may be a vector of length n.

% dipole_field returns H of size [n,3] = [Hx Hy Hz] 
% where Hx, Hy and Hz are of length n

%Demag factor d is calculated so that:
%
%  (Hx)   (dxx dxy dxz)(vMx)
%  (Hy) = (dyx dyy dyz)(vMy)
%  (Hz)   (dzx dzy dzz)(vMz)

% where vM is the dipole moment, ie M (moment per unit volume, A/m) * v (volume, m^3)

pi4i=1.0/(4.0*pi);

sep2=(x.*x+y.*y+z.*z);
sep=sqrt(sep2);
sep3=sep2.*sep;
sep5=sep3.*sep2;

sep3i=pi4i./sep3;
sep5i=3.0*pi4i./sep5;
term1x=sep5i.*x;
term1y=sep5i.*y;
term1z=sep5i.*z;

dxx=term1x.*x-sep3i;
dxy=term1x.*y;
dxz=term1x.*z;

dyx=term1y.*x;
dyy=term1y.*y-sep3i;
dyz=term1y.*z;

dzx=term1z.*x;
dzy=term1z.*y;
dzz=term1z.*z-sep3i;
