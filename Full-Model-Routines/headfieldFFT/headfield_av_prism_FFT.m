clear all
path(path,'C:\cygwin\home\kalezhi\headfield_simulations\rectangular_pole\simulations_25_july_25_by_80_by_80nm\headfields\');

hx_data = load('hx_tav_2d_data.m')*1e3; % in A/m
hy_data = load('hy_tav_2d_data.m')*1e3; % in A/m
hz_data = load('hz_tav_2d_data.m')*1e3; % in A/m
y_vec = load('y_data.m'); 

% head is aligned in the y direction: so x represents cross track
% creating shape function
% For 500 Gbits/in^2: Bx =36 nm, By = 36 nm, Lx=Bx/sqrt(2), Ly=By/sqrt(2)
% for 50 % filling factor

ao =(36/sqrt(2))/2; % half the bit length in the x direction for 500 Gbits/in^2 in nm
bo =ao; % for 2Tb/in^2 half the bit length in the y direction
t = 6;

% array filling the circle
prism_shape = shape_array_prism(ao,bo);

prism_shape
sum(sum(prism_shape))
area = 4*ao*bo


% NON UNIFORM field case
% computing the convolution
hx_av_c = hav_conv(prism_shape, hx_data, area);
hy_av_c = hav_conv(prism_shape, hy_data, area);
hz_av_c = hav_conv(prism_shape, hz_data, area);

% saving data to file to send to jim

filename_hx = 'hx_av_2d_data.m';
filename_hy = 'hy_av_2d_data.m';
filename_hz = 'hz_av_2d_data.m';
filename_y = 'y_data.m';

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

num = 151;
hx_d =hx_data(:,num); % downtrack field data along pole centre
hy_d =hy_data(:,num);
hz_d =hz_data(:,num);


% % Plotting graphs
fontsiz=40; % for intermag
linewid=10;
markersizb=22;
markersizs=16;
fontsiz_leg=20; % for intermag

vals = 1:2;
lenvals = length(vals);

Colour1=hsv2rgb([0.67*(lenvals -1)/(lenvals-1) 1 1]);
Colour2=hsv2rgb([0.67*(lenvals -2)/(lenvals-1) 1 1]);

figure(1)
clf
grid on
hold on
set(gca, 'FontSize',fontsiz);
ylabel('\it{h_{eff} (kA/m)}','FontSize',fontsiz);
xlabel('\it{down track (nm) }','FontSize',fontsiz);

heff_no =headfield_eff(hx_d, hy_d, hz_d)/1000;
heff_av_c = headfield_eff(hx_av_c(:,num), hy_av_c(:,num), hz_av_c(:,num))/1000;

plot(y_vec, heff_no,'-','LineWidth',linewid,'Color','black');
plot(y_vec, heff_av_c,'-s','LineWidth',linewid,'Color',Colour1);

hlg = legend('averaged over z only','volume averaged convolution', 2);
set(hlg,'Interpreter','tex', 'FontSize',fontsiz_leg);

figure(2)
clf
grid on
hold on
set(gca, 'FontSize',fontsiz);
ylabel('\it{h_{eff} gradient(kA/m/nm)}','FontSize',fontsiz);
xlabel('\it{down track (nm) }','FontSize',fontsiz);

heff_grad_no = gradient(heff_no,y_vec(2)-y_vec(1));
heff_av_grad_c = gradient(heff_av_c,y_vec(2)-y_vec(1));

plot(y_vec, heff_grad_no ,'-','LineWidth',linewid,'Color','black');
plot(y_vec, heff_av_grad_c,'-s','LineWidth',linewid,'Color',Colour1);

hlg = legend('averaged over z only','volume averaged convolution', 2);
set(hlg,'Interpreter','tex', 'FontSize',fontsiz_leg);

[max_grad,ind] = max(heff_av_grad_c)
opt_heff_av_c = heff_av_c(ind)*1e3 % in A/m
