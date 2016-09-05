function [havx, havy, havz] = karlqfield_av_interp(var_prop, head_prop, islandgeo_prop, h_data, r_data, s_data, s)

hg = head_prop(2);
phih = head_prop(3);

a = islandgeo_prop(2);
b = islandgeo_prop(3);

beta = b./a;
vol = islandgeo_prop(8);

rh = head_prop(7)- islandgeo_prop(6); % only absolute values considered since h1 is antisymmetric and h2 symmetric w.r.t rh

absrh = abs(rh);
h1_data = h_data(:,:,1);
h2_data = h_data(:,:,2);

if var_prop(3)==0 || var_prop(3)==3 || var_prop(3)==4 || var_prop(3)==5 % either no variations, position variations or k1 variations
    % here s_data =1, so we only perform a 1D interpolation
    h1 =  interp1(r_data, h1_data, absrh,'pchip');
    h2f = interp1(r_data, h2_data, absrh,'pchip');
else
    h1 =  interp2(r_data, s_data, h1_data, absrh, s,'spline');
    h2f = interp2(r_data, s_data, h2_data, absrh, s,'spline');
end

h1f =h1.*sign(rh); % correct signs taken into account

havx = -beta.*hg.*cos(phih).*h1f./(2*pi*vol);
havy = -beta.*hg.*sin(phih).*h1f./(2*pi*vol);
havz = beta.*hg.*h2f./(pi*vol);

end