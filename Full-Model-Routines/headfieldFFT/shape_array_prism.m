function prism_shape = shape_array_prism(ao, bo)
l_xo = 2*ao;
l_yo = 2*bo;

y_delta = l_yo - floor(l_yo);
x_delta = l_xo - floor(l_xo);

l_y = floor(l_yo) + 2*sign(y_delta);
l_x = floor(l_xo) + 2*sign(x_delta);

prism_shape = ones(l_y,l_x);

if l_x~=l_xo
    prism_shape(:,1) = x_delta/2;
    prism_shape(:,end) = x_delta/2;
end

if l_y~=l_yo
    prism_shape(1,:) = y_delta/2;
    prism_shape(end,:) = y_delta/2;
end
 
if l_y~=l_yo && l_x~=l_yo
    prism_shape(1,1) = y_delta*x_delta/4;
    prism_shape(1,end) = y_delta*x_delta/4;
    prism_shape(end,1) = y_delta*x_delta/4;
    prism_shape(end,end) = y_delta*x_delta/4;
    prism_shape(2:end-1,1) = x_delta/2;
    prism_shape(2:end-1,end) = x_delta/2;
    prism_shape(1,2:end-1) = y_delta/2;
    prism_shape(end,2:end-1) = y_delta/2;
end
end