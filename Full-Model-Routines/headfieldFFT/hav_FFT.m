function h_av = hav_FFT(shape, h_data, area)
% using FFT method
[h_lenx h_leny] = size(h_data);
shape(h_lenx, h_leny) = 0;

shape_fft = fft2(shape);
h_fft = fft2(h_data);

h_av_fft = shape_fft.*h_fft;
h_av = ifft2(h_av_fft)/area;
end