function [Hx Hy Hz] = dipole_field(Mx,My,Mz,x,y,z)

% All values in SI units
% (x,y,z) is the location of the field point relative to the dipole in m
% (Mx,My,Mz) is the dipole moment of the dipole: M.v, in (A/m)*m^3 = Am^2

% (x,y,z) may be a vector of length n.

% dipole_field returns H of size [n,3] = [Hx Hy Hz] 
% where Hx, Hy and Hz are of length n

pi4i=1.0/(4.0*pi);

sep2=(x.*x+y.*y+z.*z);
sep=sqrt(sep2);
sep3=sep2.*sep;
sep5=sep3.*sep2;

sep3i=pi4i./sep3;
term1=3.0*pi4i*(Mx.*x+My.*y+Mz.*z)./sep5;

Hx=term1.*x-Mx*sep3i;
Hy=term1.*y-My*sep3i;
Hz=term1.*z-Mz*sep3i;

