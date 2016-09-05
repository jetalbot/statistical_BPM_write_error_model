function h_extrap = extrap_h_data(zi, x_vec, y_vec, z_vec, h)
lenx = length(x_vec);
leny = length(y_vec);

stepx = x_vec(2)-x_vec(1);
stepy = y_vec(2)-y_vec(1);

h_extrap = zeros(lenx, leny);

[num_r num_c] = size(z_vec);

for i=1:stepx:lenx
    i
    for j=1:stepy:leny
        h_vec = reshape(h(j,i,:),num_r, num_c);
        h_extrap(j,i) = interp1(z_vec, h_vec, zi,'pchip','extrap');
    end
end

end