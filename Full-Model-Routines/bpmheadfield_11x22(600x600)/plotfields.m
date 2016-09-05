
%load files
% 
% hx   = load('bpm11x22_hy_av_2d_data.m');
% hy   = load('bpm11x22_hy_av_2d_data.m');
% hz   = load('bpm11x22_hz_av_2d_data.m');
% y = load('bpm11x22_y_data.m');

clear all;

hx_h_used = dlmread('bpm11x22_hx_av_2d_h_data.m'); 
hx_s_used = dlmread('bpm11x22_hx_av_2d_s_data.m'); 
hx_used   = dlmread('bpm11x22_hy_av_2d_data.m');

hy_h_used = dlmread('bpm11x22_hy_av_2d_h_data.m'); 
hy_s_used = dlmread('bpm11x22_hy_av_2d_s_data.m'); 
hy_used   = dlmread('bpm11x22_hy_av_2d_data.m');

hz_h_used = dlmread('bpm11x22_hz_av_2d_h_data.m');
hz_s_used = dlmread('bpm11x22_hz_av_2d_s_data.m');
hz_used   = dlmread('bpm11x22_hz_av_2d_data.m');
y_used = dlmread('bpm11x22_y_data.m');





%create arrays
h_eff = zeros(1,601);
h_eff_grad = zeros(1,601);
h_eff_s = zeros(1,601);
h_eff_grad_s = zeros(1,601);
h_eff_h = zeros(1,601);
h_eff_grad_h = zeros(1,601);
h_effpos = zeros(1,601);
h_effmax = zeros(1,601);



for i = 1:601

hx_h(i) = hx_h_used(i,300);
hx_s(i) = hx_s_used(i,300);
hx(i)   = hx_used(i,300);  

hy(i)   = hy_used(i,300);
hy_h(i) = hy_h_used(i,300); 
hy_s(i) = hy_s_used(i,300);

hz_h(i) = hz_h_used(i,300);
hz_s(i) = hz_s_used(i,300);
hz(i)   = hz_used(i,300);
y(i) = y_used(i)
    
    
h_eff(i) = headfield_eff(hx(i), hy(i), hz(i));
% y_spacing = 1; %y(2)-y(1); % in nm
h_eff_grad = gradient(h_eff,1);  

h_eff_s(i)= headfield_eff(hx_s(i), hy_s(i), hz_s(i));
h_eff_grad_s = gradient(h_eff_s,1);      

h_eff_h(i)= headfield_eff(hx_h(i), hy_h(i), hz_h(i));
h_eff_grad_h = gradient(h_eff_h,1); 

  hv(i) = sqrt(hx(i).^2 + hy(i).^2 + hz(i).^2);
  
    hratioz(i) = hz(i)./hv(i);
    hratioyx(i) = hy(i)./hx(i);
    theta(i) = hratioz(i)*(180/pi);
    phih(i) = hratioyx(i)*(180/pi);     
    filename = strcat('allvalues_energybarrier.m'); 
    fid=fopen(filename,'a+');
    fprintf(fid,'%g %g %g %g\n',[y(i) h_eff(i) theta(i) phih(i)]);
    fclose(fid);  
    for j = 1: 601;
    hmag(i,j) = sqrt(hx_h_used(i,j).^2 + hy_h_used(i,j).^2 + hz_h_used(i,j).^2);
    end    
end

h_eff_300 = h_eff
    h_eff_g = h_eff_grad




%theta = [79.1 78.75 77.9 77.2 75.95 74.5 72.35 70.95 68.75 66.25 63.5 60.65 57.35 53.5 49.13 42.8];
%appfield = [1.029e6 1.021e6 1.006e6 9.863e5 9.616e5 9.348e5 9.066e5 8.776e5 8.484e5 8.196e5 7.921e5 7.659e5 7.413e5 7.185e5 6.975e5 6.776e5];
% % Plotting graphs
fontsiz=20; % for intermag
linewid=5;
markersizb=15;
markersizs=10;
fontsiz_leg=10; % for intermag

figure(1)
clf
grid on
hold on
set(gca, 'FontSize',fontsiz);
ylabel('\it{h_{eff} (A/m)}','FontSize',fontsiz);
xlabel('\it{down track (nm) }','FontSize',fontsiz);
axis([0 300 -1e5 12e5])
plot(y, h_eff,'-','LineWidth',linewid,'Color','blue');
plot(y, h_eff_grad,'--','LineWidth',linewid,'Color','blue');
% plot(y, h_eff_s,'-','LineWidth',linewid,'Color','red');
% plot(y, h_eff_grad_s,'--','LineWidth',linewid,'Color','red');
% plot(y, h_eff_s,'-','LineWidth',linewid,'Color','green');
% plot(y, h_eff_grad_s,'--','LineWidth',linewid,'Color','green');


figure(2)
clf
grid on
hold on
set(gca, 'FontSize',fontsiz);
contourf(y,y,hmag)
% 
% 
% figure(3)
% clf
% grid on
% hold on
% set(gca, 'FontSize',fontsiz);
% surf(y_used,y_used,hx_used)
% 
% % 
% figure(4)
% clf
% grid on
% hold on
% set(gca, 'FontSize',fontsiz);
% surf(y,y,hv)