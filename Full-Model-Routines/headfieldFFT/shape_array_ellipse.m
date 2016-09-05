function [cyl_shape i_index j_index] = shape_array_ellipse(ao, bo)
if mod(ao,2)==0 || mod(ao,2)==1 % ao either even or odd
    if mod(bo,2)==0 || mod(bo,2)==1 % bo either even or odd
        cyl_shape = zeros(2*floor(bo)+1,2*floor(ao)+1);
    else % bo neither even or odd
        cyl_shape = zeros(2*floor(bo)+3,2*floor(ao)+1);
    end
elseif mod(bo,2)==0 || mod(bo,2)==1 % bo either even or odd
    if mod(ao,2)==0 || mod(ao,2)==1 % ao either even or odd
        cyl_shape = zeros(2*floor(bo)+1,2*floor(ao)+1);
    else % ao neither even or odd
        cyl_shape = zeros(2*floor(bo)+1,2*floor(ao)+3);
    end
else
    cyl_shape = zeros(2*floor(bo)+3,2*floor(ao)+3);
end

[ly,lx]=size(cyl_shape);

i_index = zeros(1,lx*lx);
j_index = zeros(1,ly*ly);

k=1;
a = min(ao,bo);

for i=1:lx
    for j=1:ly
        xo = i - round(lx/2);
        yo = j - round(ly/2);
        rsqrt = sqrt(((a/ao)*xo)^2 + ((a/bo)*yo)^2);
        if rsqrt + sqrt(2)/2 <= a %distance from centre to point plus radius of small circle that encloses the unit square
            cyl_shape(j,i)=1; % the unit square lies inside the ellipse
        else
            i_index(k) = i;
            j_index(k) = j;
            k=k+1;
            cyl_shape(j,i)= 0; 
        end
    end
end

i_index = i_index(i_index~=0);
j_index = j_index(j_index~=0);
end