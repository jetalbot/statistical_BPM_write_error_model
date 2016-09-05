function [nxx nyy nzz] = demagfactors_interp(x, demagdata)
% load demagdata.m
xdat = demagdata(:,1);
nxxdat = demagdata(:,2);
nyydat = demagdata(:,3);
nzzdat = demagdata(:,4);

nxx = interp1(xdat,nxxdat, x,'pchip');
nyy = interp1(xdat,nyydat, x,'pchip');
nzz = interp1(xdat,nzzdat, x,'pchip');

% pxx = [  0.00051150	 -0.00073066	 -0.01331471	  0.07375702	 -0.19840441	  0.37798433	 -0.01161559	];
% pyy = [  0.03003812	 -0.21744900	  0.67479731	 -1.18292327	  1.31590466	 -1.02385747	  0.63167716	];
% pzz = [ -0.03065445	  0.21854877	 -0.66165223	  1.10832085	 -1.11607902	  0.64502166	  0.38012085	];
% 
% nxx = polyval(pxx,x);
% nyy = polyval(pyy,x);
% nzz = polyval(pzz,x);
end