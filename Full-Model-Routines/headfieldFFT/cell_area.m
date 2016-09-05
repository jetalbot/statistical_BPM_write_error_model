function area_factor = cell_area(ro, i_index, j_index, lxo, lyo, scale_factor)
unit_cell = zeros(scale_factor);
[lx,ly]=size(unit_cell);
ro = scale_factor*ro;
area_factor = zeros(size(i_index));

xo_vec = scale_factor*(i_index - round(lxo/2));
yo_vec = scale_factor*(j_index - round(lyo/2));

for k=1:length(i_index)
    unit_cell = zeros(scale_factor);
    for i=1:lx
        for j=1:ly           
            xo = abs(xo_vec(k)) + lx/2 - i; % so that xo ranges from abs(xo_vec(k)) - lx/2 to abs(xo_vec(k)) + lx/2 -1
            yo = abs(yo_vec(k)) + ly/2 - j; % so that yo ranges from abs(yo_vec(k)) - ly/2 to abs(xo_vec(k)) + ly/2 -1
            rsqrt = sqrt(xo^2 + yo^2);
            
            if rsqrt + sqrt(2)/2 <= ro % distance from centre to point plus radius of small cirlce that encloses the unit square
                unit_cell(j,i)=1; % the unit square lies inside the circle
            else
                unit_cell(j,i)= 0;
            end
        end
    end
    area_factor(k) = sum(sum(unit_cell))/(lx*ly);
end
end