function [mean_par sigma_par] = dipolefield ( a, t, Ms)
% Calculates dipole fields for an array of island/particles distributed in a regular cubic lattice 
% in the x-y plane.

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

spacing = 30e-9; % The spacing between islands (island period or bit period). 
% The spacing is assumed to be the same in x and y (down-track and cross-track)
spacing_3_i=1.0/(spacing*spacing*spacing);

diameter = a*(10^-9); % island diameter in m 
height = t*(10^-9); % island height

Ms = Ms; % Ms in A/m

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

%All the rest is for plotting purposes only.

%Plots pretty pictures that are in some way representative.

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
mean_par = mean(hz1d_t);
disp('standard deviation of field = ')
sigma_par = std(hz1d_t);
disp(std(hz1d_t))
disp(std(hz_av))
px = normpdf(hz1d_t, mean(hz1d_t), std(hz1d_t));



% figure(4)
% hist(hz1d,size)

