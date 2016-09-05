function hdatagen(phih, gapsize, polesize, rh, zh, a, beta, t, alpha, x1, x2, npts)
xv = union(1,linspace(x1, x2, npts));
av = a./sqrt(beta*xv);
betav = beta.*xv;
tv = t.*ones(size(xv));
alphav = alpha.*ones(size(xv));

h1 = zeros(size(rh));
h2 = zeros(size(rh));
hv1 = zeros(length(xv),length(rh));
hv2 = zeros(length(xv),length(rh));
for i=1:length(xv)
    for j =1: length(rh)
        [h1(j) h2(j)] = hIntegralcone(phih, gapsize, polesize, rh(j), zh, av(i), betav(i), tv(i), alphav(i))
    end
    hv1(i,:)=h1;
    hv2(i,:)=h2;   
end

filename = 'xdata.m';
fid=fopen(filename,'w');
for i=1:length(xv)
    fprintf(fid,'%12.8f\n',xv(i));
end
fclose(fid);

filename = 'rdata.m';
fid=fopen(filename,'w');
for i=1:length(rh)
    fprintf(fid,'%12.8f\n',rh(i));
end
fclose(fid);

filename = 'h1data.m';
fid=fopen(filename,'w');
for i=1:length(xv)
    for j=1:length(rh)
        fprintf(fid,'%-12.8f\t',hv1(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);

filename = 'h2data.m';
fid=fopen(filename,'w');
for i=1:length(xv)
    for j=1:length(rh)
        fprintf(fid,'%12.8f\t',hv2(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);
end