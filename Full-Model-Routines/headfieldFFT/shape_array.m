function [cyl_shape i_index j_index] = shape_array(ro)
if mod(ro,2)==0|| mod(ro,2)==1 % either even or odd
    cyl_shape = zeros(2*floor(ro)+1,2*floor(ro)+1);
else
    cyl_shape = zeros(2*floor(ro)+3,2*floor(ro)+3);
end

[lx,ly]=size(cyl_shape);

i_index = zeros(1,lx*lx);
j_index = zeros(1,ly*ly);

k=1;
for i=1:lx
    for j=1:ly
        xo = i - round(lx/2);
        yo = j - round(ly/2);
        rsqrt = sqrt(xo^2 + yo^2);
        if rsqrt + sqrt(2)/2 <= ro % distance from centre to point plus radius of small cirlce that encloses the unit square 
            cyl_shape(j,i)=1; % the unit square lies inside the circle
        else
            i_index(k) = i;
            j_index(k) = j;
            k=k+1;
            cyl_shape(j,i)= 0; % exp(-rsqr^2);
        end
    end
end

i_index = i_index(i_index~=0);
j_index = j_index(j_index~=0);
end