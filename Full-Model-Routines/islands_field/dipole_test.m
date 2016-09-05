
%create a regular array of 20 x 20 points covering [-1:1] in x and y 
[x,y] = meshgrid(0.05:0.05:1, 0.05:0.05:1);

% z (height) = 0 everywhere
z=zeros(size(x));

errmax=0.0;
for tries=0:1000

Mx=2*(rand-0.5);
My=2*(rand-0.5);
Mz=2*(rand-0.5);

%Calculate the dipole field at positions (x,y,z) relative to a dipole 
[Hx,Hy,Hz]=dipole_field(Mx,My,Mz,x,y,z);

[Dxx Dxy Dxz Dyx Dyy Dyz Dzx Dzy Dzz]=dipole_field_d(x,y,z);

Hx2=Dxx.*Mx + Dxy.*My + Dxz.*Mz;
Hy2=Dyx.*Mx + Dyy.*My + Dyz.*Mz;
Hz2=Dzx.*Mx + Dzy.*My + Dzz.*Mz;

Habs=sqrt(Hx.*Hx+Hy.*Hy+Hz.*Hz);

Hxd=(Hx2-Hx)./Habs;
Hyd=(Hy2-Hy)./Habs;
Hzd=(Hz2-Hz)./Habs;

Herr=sqrt(Hxd.*Hxd+Hyd.*Hyd+Hzd.*Hzd);
perc_err=max(max(Herr))*100;

if perc_err > errmax
    errmax=perc_err;
end

end

errmax