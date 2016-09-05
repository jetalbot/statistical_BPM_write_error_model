function demagdatagen(a, b, t, alpha, islandgeo, varparameter, x1, x2, npts, tol)

switch varparameter
    case 1
        disp('shape variation')
        xv = union(1,linspace(x1, x2, npts));
        av = a./sqrt(xv);
        bv = b.*sqrt(xv);
        tv = t.*ones(size(xv));
    case 2
        disp('size variation')
        xv = union(1,linspace(x1, x2, npts));
        av = a.*sqrt(xv);
        bv = b.*sqrt(xv);
        tv = t.*ones(size(xv));
    otherwise
        disp('neither shape nor size')
        xv = 1;
        av = a;
        bv = b;
        tv = t;
end

switch islandgeo
    case 1
        disp('truncated elliptic cone')
        alphav = alpha.*ones(size(xv));
        [ nxx nyy nzz ] = conedemagIntegral(av, bv, tv, alphav, tol)
    case 2
        disp('elliptic cylinder')
        [ nxx nyy nzz ] = cyldemagIntegral(av, bv, tv, tol)
    case 3
        disp('prism')
        [nxx nyy nzz] = prismdemagIntegral(av, bv, tv/2) % since c =t/2
    otherwise
        disp('Unknown geometry.')
end

filename = 'demagdata.m';
fid=fopen(filename,'w');

for i=1:length(xv)
    fprintf(fid,'%12.8f %12.8f %12.8f %12.8f\n',[xv(i); nxx(i); nyy(i); nzz(i)]);
end
fclose(fid);
end