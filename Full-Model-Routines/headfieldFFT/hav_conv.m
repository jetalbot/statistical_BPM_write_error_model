function h_av = hav_conv(shape, h_data, area)
% direct convolution
h_av = conv2(h_data, shape, 'same')/area;
end