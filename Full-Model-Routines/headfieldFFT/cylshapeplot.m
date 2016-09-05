% creating shape function
ro =6.7/2; 
t = 10;
% array filling the circle
[cyl_shape i_index j_index] = shape_array(ro);

cyl_shape 
sum(sum(cyl_shape))

% magnifying each unit square to get the area fraction

[lx,ly]=size(cyl_shape);
k_vec = 1:10:500

area_vec = zeros(size(k_vec));
cell_magn = zeros(size(k_vec));
area = pi*ro^2;

scale_factor =0;
for k=1:50
    k
    scale_factor = scale_factor + 10
    cell_magn(k) = scale_factor;
    area_factor = cell_area(ro, i_index, j_index, lx, ly, scale_factor);
    for i=1:length(area_factor)
        cyl_shape(j_index(i),i_index(i)) = area_factor(i);
    end
    % cyl_shape
    area_vec(k) = sum(sum(cyl_shape));

end
% 
% % Plotting graphs
fontsiz=40; % for intermag
linewid=10;
markersizb=22;
markersizs=16;
fontsiz_leg=35; % for intermag

figure(1)
clf
grid on
hold on
set(gca, 'FontSize',fontsiz);
ylabel('\it{area}','FontSize',fontsiz);
xlabel('\it{cell magnification}','FontSize',fontsiz);
plot(cell_magn, area*ones(size(cell_magn)),'-','LineWidth',linewid,'Color','black'); 
plot(cell_magn, area_vec,'-','LineWidth',linewid,'Color','blue'); 

% hlg = legend('volume averaged direct','area averaged conv','area averaged FFT','no averaging', 3);
% set(hlg,'Interpreter','tex', 'FontSize',fontsiz_leg);
