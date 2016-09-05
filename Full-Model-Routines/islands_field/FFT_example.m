
% Calculates dipole fields for an array of island/particles distributed in a regular cubic lattice 
% in the x-y plane.

fontsiz=18;
linewid=3.0;
hollow_marker=12;
filled_marker=24;

% Uses SI units throughout

% Only calculates Hz(Mz) (assumes perpendicular anisotropy and that because
% there is no z-height modulation the Hx, Hy components all average out
% to zero.

% Since there is no z-height difference the position of one island relative
% to another is just (x,y) = (delta_i,delta_j)*spacing, r=sqrt(x^2+y^2) and Hz = Mz/((4.0*pi)*r^3)

% Define the size of 2-D arrays
size=1024;
% middle of the arrays is:
s2=size/2 + 1;

spacing = 18e-9; % The spacing between islands (island period or bit period). 
% The spacing is assumed to be the same in x and y (down-track and cross-track)
spacing_3_i=1.0/(spacing*spacing*spacing);

diameter = 8e-9; % island diameter
height = 6e-9; % island height

Ms = 1000e3 % Ms in A/m

Moment = Ms*(pi*diameter^2/4)*height; %note that if you want to get clever you can cancel the pi here with the one in d(ij)

% Set up the m array to be randomly +1 or -1 (ie Mz is either up or down):
mz = Moment * 2 * ((rand(size,size)>0.5) - 0.5);

%Set up the d array
dzz = zeros(size,size);
for i=1:size
    for j=1:size
        if i==s2 & j==s2
            dzz(i,j)=0; %NOTE - self-interaction is excluded. 
                      %It's assumed that the internal (self) demag field is dealt with in shape anisotropy. 
        else
            dzz(i,j)=-1./((4*pi)*((i-s2)*(i-s2)+(j-s2)*(j-s2))^(3/2));
        end
    end
end

dzz = dzz * spacing_3_i;

%Calculate h at every island due to all other islands using the FFT method
Mz=fft2(fftshift(mz));
Dzz=fft2(fftshift(dzz));
Hz=Mz.*Dzz;
hz=fftshift(ifft2(Hz)); %in A/m

%Done

% figure(1)
% pcolor(mz)
% colorbar
% 
% % Plot abs(dzz) on a log scale:
% figure(2)
% dplot=log10(abs(dzz));
% pcolor(dplot)
% colorbar
% 
% figure(3)
% pcolor(hz)
% colorbar

%Print out the statistics
hz1d=reshape(hz,size*size,1);
hz1d_t = sort(hz1d);
hz_av = (hz1d_t-mean(hz1d_t))/std(hz1d_t);
disp('mean field = ')
disp(mean(hz1d_t))
disp(mean(hz_av))
disp('standard deviation of field = ')
disp(std(hz_av))
disp(std(hz1d_t))


figure(4)
[nh,xh]=hist(hz1d_t,50);
hold on
plot(xh,nh/(sum(nh)*(xh(2)-xh(1))),'b');
plot(hz1d_t, normpdf(hz1d_t, mean(hz1d_t), std(hz1d_t)),'r','LineWidth',linewid,'MarkerSize',hollow_marker)


x=-10:0.1:10
figure(5)
hold on
plot(x, normpdf(x, (mean(hz_av)), std(hz_av)),'r','LineWidth',linewid,'MarkerSize',hollow_marker)
plot(x, normpdf(x, 0, 0.075),'--b','LineWidth',linewid,'MarkerSize',hollow_marker)

%plot((hz1d-mean(hz1d))/std(hz1d),normpdf(hz1d, mean(hz1d), std(hz1d)),'b','LineWidth',linewid,'MarkerSize',hollow_marker)
