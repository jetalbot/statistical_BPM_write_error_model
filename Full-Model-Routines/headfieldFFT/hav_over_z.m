function h_av_oz = hav_over_z(shape, h_data, lx, ly, t)
h_av_oz = zeros(lx,ly);

[row,col,shape_nz]=find(shape); % find non zero elements of shape and store them in shape_nz
len_shape = length(shape_nz);

for i=1:lx
    for j=1:ly
        h = reshape(h_data(i,j,:),1,[]);
        % h_conv = conv(h, shape); % alternatively can use conv(h, shape_nz) and get the same answer
        h_conv = conv(h(col), shape_nz); % used h corresponding to shape_nz since there are leading zeros in shape
        h_av_oz(i,j) = h_conv(len_shape)/t; % since the h_conv(len_shape) gives the right result
    end
end
end