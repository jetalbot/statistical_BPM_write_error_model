%% HEAD FIELD VALUES INTERPOLATION
% Extracts the volume averaged headfield files provided (for volume
% averaging see head_gen.m) and interpolates for the current head position.
% Provides headfield values for all axes, for hard and soft layers or for
% single layers. 
%
%%
function [hx_av hy_av hz_av hx_av_h hy_av_h hz_av_h hx_av_s hy_av_s hz_av_s] = real_headfield_av_vec_interp(var_prop, head_prop, islandgeo_prop, h_data, x_data, y_data, s_data, s, h_data_h, h_data_s)

hg = head_prop(2); % head gap field 
y = head_prop(11) + islandgeo_prop(6) - head_prop(7); % down track position
x = head_prop(12) + islandgeo_prop(7) - head_prop(8); % cross track position

% The head field are assumed to already be averaged over island volume, so
% there are 'slices' of field values for several values of parameter z, to
% interpolate in 3D. 
            
if var_prop(3)==0 || var_prop(3)==2 ||var_prop(3)==3 || var_prop(3)==4 || var_prop(3)==5
    % either no variations, position variations or k1 variations
    
    % if length(islandgeo_prop)>10 % here we have ECC structure
    hx_data = h_data(:,:,1);
    hy_data = h_data(:,:,2);
    hz_data = h_data(:,:,3);
    
    hx_data_h = h_data_h(:,:,1);
    hy_data_h = h_data_h(:,:,2);
    hz_data_h = h_data_h(:,:,3);
    
    hx_data_s = h_data_s(:,:,1);
    hy_data_s = h_data_s(:,:,2);
    hz_data_s = h_data_s(:,:,3);
    
    if length(x_data)==1
        % If x is one, there is no cross track field to consider, only 1
        % scan line, so we perform a 1D interpolation
        hx_av =  interp1(y_data, hx_data, y,'pchip'); % verify
        hy_av =  interp1(y_data, hy_data, y,'pchip');
        hz_av =  interp1(y_data, hz_data, y,'pchip');
        
        hx_av_h =  interp1(y_data, hx_data_h, y,'pchip'); % verify
        hy_av_h =  interp1(y_data, hy_data_h, y,'pchip');
        hz_av_h =  interp1(y_data, hz_data_h, y,'pchip');
        
        hx_av_s =  interp1(y_data, hx_data_s, y,'pchip'); % verify
        hy_av_s =  interp1(y_data, hy_data_s, y,'pchip');
        hz_av_s =  interp1(y_data, hz_data_s, y,'pchip');
        
    else
        % here length(x_data) > 1, z_data =1, so we only perform a 2D
        % interpolation
        hx_av =  interp2(x_data, y_data, hx_data, x, y,'pchip'); % verify
        hy_av =  interp2(x_data, y_data, hy_data, x, y,'pchip');
        hz_av =  interp2(x_data, y_data, hz_data, x, y,'pchip');
        
        hx_av_h =  interp2(x_data, y_data, hx_data_h, x, y,'pchip'); % verify
        hy_av_h =  interp2(x_data, y_data, hy_data_h, x, y,'pchip');
        hz_av_h =  interp2(x_data, y_data, hz_data_h, x, y,'pchip');
        
        hx_av_s =  interp2(x_data, y_data, hx_data_s, x, y,'pchip'); % verify
        hy_av_s =  interp2(x_data, y_data, hy_data_s, x, y,'pchip');
        hz_av_s =  interp2(x_data, y_data, hz_data_s, x, y,'pchip');
    end
else
    %Perform a 3D interpolation
    hx_data = h_data(:,:,:,1);
    hy_data = h_data(:,:,:,2);
    hz_data = h_data(:,:,:,3);
    
    hx_data_h = h_data_h(:,:,:,1);
    hy_data_h = h_data_h(:,:,:,2);
    hz_data_h = h_data_h(:,:,:,3);
    
    hx_data_s = h_data_s(:,:,:,1);
    hy_data_s = h_data_s(:,:,:,2);
    hz_data_s = h_data_s(:,:,:,3);

    hx_av =  interp3(x_data, y_data, s_data, hx_data, x, y, s,'cubic');
    hy_av =  interp3(x_data, y_data, s_data, hy_data, x, y, s,'cubic');
    hz_av =  interp3(x_data, y_data, s_data, hz_data, x, y, s,'cubic');
    
    hx_av_h =  interp3(x_data, y_data, s_data, hx_data_h, x, y, s,'cubic');
    hy_av_h =  interp3(x_data, y_data, s_data, hy_data_h, x, y, s,'cubic');
    hz_av_h =  interp3(x_data, y_data, s_data, hz_data_h, x, y, s,'cubic');
    
    hx_av_s =  interp3(x_data, y_data, s_data, hx_data_s, x, y, s,'cubic');
    hy_av_s =  interp3(x_data, y_data, s_data, hy_data_s, x, y, s,'cubic');
    hz_av_s =  interp3(x_data, y_data, s_data, hz_data_s, x, y, s,'cubic');
end

hx_av =hg.*hx_av;
hy_av =hg.*hy_av;
hz_av =hg.*hz_av;

hx_av_h =hg.*hx_av_h;
hy_av_h =hg.*hy_av_h;
hz_av_h =hg.*hz_av_h;

hx_av_s =hg.*hx_av_s;
hy_av_s =hg.*hy_av_s;
hz_av_s =hg.*hz_av_s;

end
