function [hx_avdata hy_avdata hz_avdata]=  headfield_interp2(y_data, hx_av_2d_data, hy_av_2d_data, hz_av_2d_data,head_pos)
hx_avdata =  interp2(y_data, y_data, hx_av_2d_data, max(y_data)/2+head_pos, y_data,'spline')';
hy_avdata =  interp2(y_data, y_data, hy_av_2d_data, max(y_data)/2+head_pos, y_data,'spline')';
hz_avdata =  interp2(y_data, y_data, hz_av_2d_data, max(y_data)/2+head_pos, y_data,'spline')';
end