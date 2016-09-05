function h = h_int(z, h_vec, z_vec)
h = zeros(size(z));
for i=1:length(z)
    h(i) = interp1(z_vec, h_vec, z(i),'pchip','extrap');
end
end