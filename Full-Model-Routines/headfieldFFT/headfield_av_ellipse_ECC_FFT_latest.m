clear;
path(path,'C:\Documents and Settings\kalezhi\My Documents\PhDPlus\model_routines\write_simulations\bpmheadfield_11x22nm');

% % % bpm_11_by_22 pole
% % Load the data and extract the (x,y,z, hx, hy, hz) information:
headfield = load('bpm-4tb-11x22-ss7-ts6-air5-rl6-il1-1nm.hed');

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
ao =8/2; % for 2Tb/in^2 in the down-track direction
bo =ao; % for 2Tb/in^2 in the cros-track direction
t = 6;

% lower hard layer
t_h = 4; %3.23; %5;
% begin with vertical component first, average head field over z so we end up with a 2D matrix
thicknes_shape_h = zeros(size(z_vec)); % z_vec ranges from 1 to 7 after inserting a new layer
% island z begins at 1.5 to 1.5+6
x_start = 0.5;
x_rem = t_h;
for i=1:length(thicknes_shape_h)
    if i == 1
        thicknes_shape_h(i) = x_start;
        x_rem = x_rem - x_start;
    elseif x_rem >=1
        thicknes_shape_h(i) = 1;
        x_rem = x_rem - 1;
    elseif x_rem >=0 && x_rem < 1
        thicknes_shape_h(i) = x_rem;
        x_inter_h = x_rem;
        x_rem = x_rem - 1;
        i_inter = i;
    else
        thicknes_shape_h(i) = 0;
    end
end
sum(thicknes_shape_h) % should be t_h

% upper soft layer
t_s = 2; %t - t_h; %6.77; %5;
% begin with vertical component first, average head field over z so we end up with a 2D matrix
thicknes_shape_s = zeros(size(z_vec)); % z_vec ranges from 1 to 13 after inserting a new layer
% island z begins at 1.5 to 1.5+10
x_start =  1 - x_inter_h;
x_rem = t_s;
for i=1:length(thicknes_shape_s)
    if i < i_inter
        thicknes_shape_s(i) = 0;
    elseif i == i_inter
        thicknes_shape_s(i) = x_start;
        x_rem = x_rem - x_start;
    elseif i > i_inter && x_rem >=1
        thicknes_shape_s(i) = 1;
        x_rem = x_rem - 1;
    elseif i > i_inter && x_rem >=0 && x_rem < 1
        thicknes_shape_s(i) = x_rem;
        x_rem = x_rem - 1;
    else
        thicknes_shape_s(i) = 0;
    end
end
sum(thicknes_shape_s) % should be t_s

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
hx_av_s = hav_over_z(thicknes_shape_s, hx_data, lx, ly, t_s);
hy_av_s = hav_over_z(thicknes_shape_s, hy_data, lx, ly, t_s);
hz_av_s = hav_over_z(thicknes_shape_s, hz_data, lx, ly, t_s);

hx_av_h = hav_over_z(thicknes_shape_h, hx_data, lx, ly, t_h);
hy_av_h = hav_over_z(thicknes_shape_h, hy_data, lx, ly, t_h);
hz_av_h = hav_over_z(thicknes_shape_h, hz_data, lx, ly, t_h);

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
hx_av_c_s = hav_conv(cyl_shape, hx_av_s, area);
hy_av_c_s = hav_conv(cyl_shape, hy_av_s, area);
hz_av_c_s = hav_conv(cyl_shape, hz_av_s, area);

hx_av_c_h = hav_conv(cyl_shape, hx_av_h, area);
hy_av_c_h = hav_conv(cyl_shape, hy_av_h, area);
hz_av_c_h = hav_conv(cyl_shape, hz_av_h, area);

hx_av_c = hav_conv(cyl_shape, hx_av, area);
hy_av_c = hav_conv(cyl_shape, hy_av, area);
hz_av_c = hav_conv(cyl_shape, hz_av, area);

toc

% saving data to file to send to jim

filename_hx_h = 'bpm11x22_hx_av_2d_h_data.m';
filename_hy_h = 'bpm11x22_hy_av_2d_h_data.m';
filename_hz_h = 'bpm11x22_hz_av_2d_h_data.m';

filename_hx_s = 'bpm11x22_hx_av_2d_s_data.m';
filename_hy_s = 'bpm11x22_hy_av_2d_s_data.m';
filename_hz_s = 'bpm11x22_hz_av_2d_s_data.m';

filename_hx = 'bpm11x22_hx_av_2d_data.m';
filename_hy = 'bpm11x22_hy_av_2d_data.m';
filename_hz = 'bpm11x22_hz_av_2d_data.m';
filename_y = 'bpm11x22_y_data.m';

fid_x_h =fopen(filename_hx_h,'w');
fid_y_h =fopen(filename_hy_h,'w');
fid_z_h =fopen(filename_hz_h,'w');

fid_x_s =fopen(filename_hx_s,'w');
fid_y_s =fopen(filename_hy_s,'w');
fid_z_s =fopen(filename_hz_s,'w');

fid_x =fopen(filename_hx,'w');
fid_y =fopen(filename_hy,'w');
fid_z =fopen(filename_hz,'w');
fid =fopen(filename_y,'w');

[h_len_x,h_len_y] =size(hx_av_c);

for i=1:h_len_x

    fprintf(fid,'%12.8f\n',y_vec(i));

    for j=1:h_len_y
        fprintf(fid_x_s,'%12.8f\t',hx_av_c_s(i,j));
        fprintf(fid_y_s,'%12.8f\t',hy_av_c_s(i,j));
        fprintf(fid_z_s,'%12.8f\t',hz_av_c_s(i,j));

        fprintf(fid_x_h,'%12.8f\t',hx_av_c_h(i,j));
        fprintf(fid_y_h,'%12.8f\t',hy_av_c_h(i,j));
        fprintf(fid_z_h,'%12.8f\t',hz_av_c_h(i,j));

        fprintf(fid_x,'%12.8f\t',hx_av_c(i,j));
        fprintf(fid_y,'%12.8f\t',hy_av_c(i,j));
        fprintf(fid_z,'%12.8f\t',hz_av_c(i,j));
    end
    fprintf(fid_x_s,'\n');
    fprintf(fid_y_s,'\n');
    fprintf(fid_z_s,'\n');
    
    fprintf(fid_x_h,'\n');
    fprintf(fid_y_h,'\n');
    fprintf(fid_z_h,'\n');
    
    fprintf(fid_x,'\n');
    fprintf(fid_y,'\n');
    fprintf(fid_z,'\n');
end
fclose(fid);

fclose(fid_x_s);
fclose(fid_y_s);
fclose(fid_z_s);

fclose(fid_x_h);
fclose(fid_y_h);
fclose(fid_z_h);

fclose(fid_x);
fclose(fid_y);
fclose(fid_z);

col_mid = 301;
% unprocessed data
hx_d =hx_data(:,col_mid,3); % downtrack field data for 3rd layer
hy_d =hy_data(:,col_mid,3);
hz_d =hz_data(:,col_mid,3);

% % Plotting graphs
fontsiz=40; % for intermag
linewid=10;
markersizb=22;
markersizs=16;
fontsiz_leg=20; % for intermag

vals = 1:3;
lenvals = length(vals);

Colour1=hsv2rgb([0.67*(lenvals -3)/(lenvals-1) 1 1]);
Colour2=hsv2rgb([0.67*(lenvals -2)/(lenvals-1) 1 1]);
Colour3=hsv2rgb([0.67*(lenvals -1)/(lenvals-1) 1 1]);

% figure(1)
% clf
% grid on
% hold on
% set(gca, 'FontSize',fontsiz);
% ylabel('\it{h_{eff} (A/m)}','FontSize',fontsiz);
% xlabel('\it{down track (nm) }','FontSize',fontsiz);
% 
% heff_no =headfield_eff(hx_d, hy_d, hz_d);
% 
% heff_av_z_s = headfield_eff(hx_av_s(:,col_mid), hy_av_s(:,col_mid), hz_av_s(:,col_mid));
% heff_av_z_h = headfield_eff(hx_av_h(:,col_mid), hy_av_h(:,col_mid), hz_av_h(:,col_mid));
% heff_av_z = headfield_eff(hx_av(:,col_mid), hy_av(:,col_mid), hz_av(:,col_mid));
% 
% heff_av_c_s = headfield_eff(hx_av_c_s(:,col_mid), hy_av_c_s(:,col_mid), hz_av_c_s(:,col_mid));
% heff_av_c_h = headfield_eff(hx_av_c_h(:,col_mid), hy_av_c_h(:,col_mid), hz_av_c_h(:,col_mid));
% heff_av_c = headfield_eff(hx_av_c(:,col_mid), hy_av_c(:,col_mid), hz_av_c(:,col_mid));
% % heff_av_d = headfield_eff(hx_av_data, hy_av_data, hz_av_data);
% 
% plot(y_vec, heff_no,'-','LineWidth',linewid,'Color','black');
% 
% plot(y_vec, heff_av_z,'s','LineWidth',linewid,'Color',Colour1);
% plot(y_vec, heff_av_z_s,'s','LineWidth',linewid,'Color',Colour2);
% plot(y_vec, heff_av_z_h,'s','LineWidth',linewid,'Color',Colour3);
% 
% plot(y_vec, heff_av_c,'-s','LineWidth',linewid,'Color',Colour1);
% plot(y_vec, heff_av_c_s,'-s','LineWidth',linewid,'Color',Colour2);
% plot(y_vec, heff_av_c_h,'-s','LineWidth',linewid,'Color',Colour3);
% 
% % plot(y_data, heff_av_d,'-','LineWidth',linewid,'Color',Colour3);
% 
% % hlg = legend('no averaging','averaged over z only','volume averaged convolution', 'volume averaged direct', 4);
% hlg = legend('no averaging','averaged convolution over z','averaged over z only soft','averaged over z only hard','volume averaged convolution','volume averaged convolution soft', 'volume averaged convolution hard',3);
% set(hlg,'Interpreter','tex', 'FontSize',fontsiz_leg);


figure(2)
clf
grid on
hold on
set(gca, 'FontSize',fontsiz);
ylabel('\it{h_{z} (kA/m)}','FontSize',fontsiz);
xlabel('\it{down track (nm) }','FontSize',fontsiz);

hz_no = hz_d/1e3;

hz_av_z_s_1 = hz_av_s(:,col_mid)/1e3;
hz_av_z_h_1 = hz_av_h(:,col_mid)/1e3;
hz_av_z_1 = hz_av(:,col_mid)/1e3;

hz_av_c_s_1 = hz_av_c_s(:,col_mid)/1e3;
hz_av_c_h_1 = hz_av_c_h(:,col_mid)/1e3;
hz_av_c_1 = hz_av_c(:,col_mid)/1e3;
% heff_av_d = headfield_eff(hx_av_data, hy_av_data, hz_av_data);

plot(y_vec, hz_no,'-','LineWidth',linewid,'Color','black');

plot(y_vec, hz_av_z_1,'s','LineWidth',linewid,'Color',Colour1);
plot(y_vec, hz_av_z_s_1,'s','LineWidth',linewid,'Color',Colour2);
plot(y_vec, hz_av_z_h_1,'s','LineWidth',linewid,'Color',Colour3);

plot(y_vec, hz_av_c_1,'-s','LineWidth',linewid,'Color',Colour1);
plot(y_vec, hz_av_c_s_1,'-s','LineWidth',linewid,'Color',Colour2);
plot(y_vec, hz_av_c_h_1,'-s','LineWidth',linewid,'Color',Colour3);

% plot(y_data, heff_av_d,'-','LineWidth',linewid,'Color',Colour3);

% hlg = legend('no averaging','averaged over z only','volume averaged convolution', 'volume averaged direct', 4);
hlg = legend('no averaging','averaged convolution over z','averaged over z only soft','averaged over z only hard','volume averaged convolution','volume averaged convolution soft', 'volume averaged convolution hard',3);
set(hlg,'Interpreter','tex', 'FontSize',fontsiz_leg);

figure(3)
clf
grid on
hold on
set(gca, 'FontSize',fontsiz);
ylabel('\it{h_{z} (kA/m)}','FontSize',fontsiz);
xlabel('\it{down track (nm) }','FontSize',fontsiz);

hz_no = hz_d/1e3;

hz_av_z_s_1 = hz_av_s(:,col_mid)/1e3;
hz_av_z_h_1 = hz_av_h(:,col_mid)/1e3;
hz_av_z_1 = hz_av(:,col_mid)/1e3;

hz_av_c_s_1 = hz_av_c_s(:,col_mid)/1e3;
hz_av_c_h_1 = hz_av_c_h(:,col_mid)/1e3;
hz_av_c_1 = hz_av_c(:,col_mid)/1e3;
% heff_av_d = headfield_eff(hx_av_data, hy_av_data, hz_av_data);

plot(y_vec, -hz_no,'-','LineWidth',linewid,'Color','black');

plot(y_vec, -hz_av_z_1,'s','LineWidth',linewid,'Color',Colour1);
plot(y_vec, -hz_av_z_s_1,'s','LineWidth',linewid,'Color',Colour2);
plot(y_vec, -hz_av_z_h_1,'s','LineWidth',linewid,'Color',Colour3);

plot(y_vec, -hz_av_c_1,'-s','LineWidth',linewid,'Color',Colour1);
plot(y_vec, -hz_av_c_s_1,'-s','LineWidth',linewid,'Color',Colour2);
plot(y_vec, -hz_av_c_h_1,'-s','LineWidth',linewid,'Color',Colour3);

% plot(y_data, heff_av_d,'-','LineWidth',linewid,'Color',Colour3);

% hlg = legend('no averaging','averaged over z only','volume averaged convolution', 'volume averaged direct', 4);
hlg = legend('no averaging','averaged convolution over z','averaged over z only soft','averaged over z only hard','volume averaged convolution','volume averaged convolution soft', 'volume averaged convolution hard',3);
set(hlg,'Interpreter','tex', 'FontSize',fontsiz_leg);

hz_no_grad = gradient(-hz_no);
hz_av_z_s_1_grad = gradient(-hz_av_s(:,col_mid)/1e3);
hz_av_z_h_1_grad = gradient(-hz_av_h(:,col_mid)/1e3);
hz_av_z_1_grad = gradient(-hz_av(:,col_mid)/1e3);

hz_av_c_s_1_grad = gradient(-hz_av_c_s(:,col_mid)/1e3);
hz_av_c_h_1_grad = gradient(-hz_av_c_h(:,col_mid)/1e3);
hz_av_c_1_grad = gradient(-hz_av_c(:,col_mid)/1e3);

figure(4)
clf
grid on
hold on
set(gca, 'FontSize',fontsiz);
ylabel('\it{grad (h_{z} (kA/m)/nm}','FontSize',fontsiz);
xlabel('\it{down track (nm) }','FontSize',fontsiz);

plot(y_vec, hz_no_grad,'-','LineWidth',linewid,'Color','black');

plot(y_vec, hz_av_z_s_1_grad,'s','LineWidth',linewid,'Color',Colour1);
plot(y_vec, hz_av_z_h_1_grad,'s','LineWidth',linewid,'Color',Colour2);
plot(y_vec, hz_av_z_1_grad,'s','LineWidth',linewid,'Color',Colour3);

plot(y_vec, hz_av_c_s_1_grad,'-s','LineWidth',linewid,'Color',Colour1);
plot(y_vec, hz_av_c_h_1_grad,'-s','LineWidth',linewid,'Color',Colour2);
plot(y_vec, hz_av_c_1_grad,'-s','LineWidth',linewid,'Color',Colour3);

% plot(y_data, heff_av_d,'-','LineWidth',linewid,'Color',Colour3);

% hlg = legend('no averaging','averaged over z only','volume averaged convolution', 'volume averaged direct', 4);
hlg = legend('no averaging','averaged convolution over z','averaged over z only soft','averaged over z only hard','volume averaged convolution','volume averaged convolution soft', 'volume averaged convolution hard',3);
set(hlg,'Interpreter','tex', 'FontSize',fontsiz_leg);

