clear;
path(path,'C:\Documents and Settings\kalezhi\My Documents\PhDPlus\model_routines\write_simulations\bpm_13_by_40');

% % % bpm_11_by_22 pole
% % Load the data and extract the (x,y,z, hx, hy, hz) information:
headfield = load('bpm-13x40-th8-ss10-ts5-1nm.hed');

% get required reshaped head field data
% the function also adds a missing layer, if this is not necessary, then
% set z_extrap to -1;

z_extrap = 1; % add a new layer at z = 1
[hx_data hy_data hz_data x_vec y_vec z_vec] = head_gen_latest(headfield, z_extrap);
clear headfield

% the fields are usually in Oesterds, to convert field to SI units, A/m we multiply the fields by 1e3/(4*pi).

to_apm =1e3/(4*pi); % since 1Oe = 10^3/(4*pi) A/m

hx_data = hx_data*to_apm;
hy_data = hy_data*to_apm;
hz_data = hz_data*to_apm;

% creating shape function
ao =12.5/2; % semi-minor/minor axis in the down-track direction
bo =ao; % semi-major/minor axis in the cros-track direction
t = 10;

% begin with vertical component first, average head field over z so we end up with a 2D matrix
thicknes_shape = zeros(size(z_vec)); % z_vec ranges from 1 to 13 after inserting a new layer
% island z begins at 1.5 to 1.5+10
for i=1:length(thicknes_shape)
    if i == 1
        thicknes_shape(i) = 0.5;
    elseif i >1 && i < round(t)+1
        thicknes_shape(i) = 1;
    elseif i==round(t)+1
        thicknes_shape(i) = 0.5;
    else
        thicknes_shape(i) = 0;
    end
end
% sum(thicknes_shape) % should be t

% averaging over z first
lx = length(x_vec);
ly = length(y_vec);
tic

hx_av = hav_over_z(thicknes_shape, hx_data, lx, ly, t);
hy_av = hav_over_z(thicknes_shape, hy_data, lx, ly, t);
hz_av = hav_over_z(thicknes_shape, hz_data, lx, ly, t);

% array filling the circle
[cyl_shape i_index j_index] = shape_array_ellipse(ao,bo);

cyl_shape
sum(sum(cyl_shape))
area = pi*ao*bo

% magnifying each unit square to get the area fraction
scale_factor = 500;
[ly,lx]=size(cyl_shape);
area_factor = cell_area_ellipse(ao, bo, i_index, j_index, lx, ly, scale_factor)

% updating the flagged unit squares
for i=1:length(area_factor)
    cyl_shape(j_index(i),i_index(i)) = area_factor(i);
end
cyl_shape
sum(sum(cyl_shape))
area = pi*ao*bo

% NON UNIFORM field case
% computing the convolution

hx_av_c = hav_conv(cyl_shape, hx_av, area);
hy_av_c = hav_conv(cyl_shape, hy_av, area);
hz_av_c = hav_conv(cyl_shape, hz_av, area);

toc

% saving data to file

% saving data to file
dir_name ='bpm_headfield_averaged_files'; % directory name

mkdir(dir_name); % create directory

cd(dir_name);

filename_hx = 'bpm_headfield_hx_av_2d_data.m';
filename_hy = 'bpm_headfield_hy_av_2d_data.m';
filename_hz = 'bpm_headfield_hz_av_2d_data.m';
filename_y = 'bpm_headfield_y_data.m';

fid_x =fopen(filename_hx,'w');
fid_y =fopen(filename_hy,'w');
fid_z =fopen(filename_hz,'w');
fid =fopen(filename_y,'w');

[h_len_x,h_len_y] =size(hx_av_c);

for i=1:h_len_x

    fprintf(fid,'%12.8f\n',y_vec(i));

    for j=1:h_len_y
        fprintf(fid_x,'%12.8f\t',hx_av_c(i,j));
        fprintf(fid_y,'%12.8f\t',hy_av_c(i,j));
        fprintf(fid_z,'%12.8f\t',hz_av_c(i,j));
    end
    
    fprintf(fid_x,'\n');
    fprintf(fid_y,'\n');
    fprintf(fid_z,'\n');
end
fclose(fid);

fclose(fid_x);
fclose(fid_y);
fclose(fid_z);

cd('..');