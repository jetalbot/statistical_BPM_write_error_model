function [hx_avdata hy_avdata  hz_avdata ] = headfield_interp(x_data, y_data, hx_av_data, hy_av_data, hz_av_data, head_centre, head_pos)

[s1 s2] = size(hx_av_data);

if min(s1,s2)==1 % here we have
    hx_avdata =  interp1(y_data, hx_av_data, y_data,'spline')';
    hy_avdata =  interp1(y_data, hy_av_data, y_data,'spline')';
    hz_avdata =  interp1(y_data, hz_av_data, y_data,'spline')';
    
elseif min(s1,s2) >1 % here we have a 2-D array
    hx_avdata =  interp2(x_data, y_data, hx_av_data,  head_centre+head_pos, y_data,'spline')';
    hy_avdata =  interp2(x_data, y_data, hy_av_data,  head_centre+head_pos, y_data,'spline')';
    hz_avdata =  interp2(x_data, y_data, hz_av_data,  head_centre+head_pos, y_data,'spline')';
else
    disp('invalid head field distribution.')
end
end