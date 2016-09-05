 function headfielddatagen(realhead_pos_prop, realhead_field_prop, head_prop, islandgeo_prop, tperiod, tptiny_prop, var_prop, tol_prop, step_vect, s1, s2, npts, num_steps)
a = islandgeo_prop(2);
b = islandgeo_prop(3);
t = islandgeo_prop(4);
alpha = islandgeo_prop(5);

tol = tol_prop(3);

switch var_prop(3) % varparameter
    case 1
        disp('shape variation')
        sv = union(1,linspace(s1, s2, npts));
        av = a./sqrt(sv);
        bv = b.*sqrt(sv);
        tv = t.*ones(size(sv));
        alphav = alpha.*ones(size(sv));
    case 2
        disp('size variation')
        sv = union(1,linspace(s1, s2, npts));
        av = a.*sqrt(sv);
        bv = b.*sqrt(sv);
        tv = t.*ones(size(sv));
        alphav = alpha.*ones(size(sv));
    otherwise
        disp('neither shape nor size; this works either for no variations, position variations, k1 variations, ms variations')
        sv = 1;
        av = a;
        bv = b;
        tv = t;
        alphav = alpha;
end


switch head_prop(1) % head type
    case 1 % headtype = 1 represents karlqvist type
        disp('headtype = 1: karlqvist type')
        
        tptiny = linspace(0,tptiny_prop(1),tptiny_prop(2));
        vel = head_prop(9);
        rh = vel*tperiod.*tptiny; % separation between island and head

        phih = head_prop(3);
        gapsize = head_prop(4);
        polesize = head_prop(5);
        zh = head_prop(6);

        betav = bv./av;

        h1 = zeros(length(sv),length(rh));
        h2 = zeros(length(sv),length(rh));
        lens =length(sv);
        switch islandgeo_prop(1)
            case 1
                disp('truncated elliptic cone')

                for i=1:lens
                    for j =1: length(rh)
                        [h1(i,j) h2(i,j)] = hIntegralcone(phih, gapsize, polesize, rh(j), zh, av(i), betav(i), tv(i), alphav(i), tol)
                    end
                end
            case 2
                disp('elliptic cylinder')
                for i=1:lens
                    for j =1: length(rh)
                        [h1(i,j) h2(i,j)] = hIntegralcyl(phih, gapsize, polesize, rh(j), zh, av(i), betav(i), tv(i), tol)
                    end
                end
            case 3
                disp('prism')
                dv = gapsize/2 -zh -tv;
                for i=1:lens
                    for j =1: length(rh)
                        [h1(i,j) h2(i,j)] = hIntegralpris(phih, gapsize, polesize, dv(i), rh(j), av(i), bv(i), tv(i), tol)
                    end
                end
            otherwise
                disp('Unknown geometry.')
        end

        filename = 'h1data.m';
        fid=fopen(filename,'w');
        for i=1:lens
            for j=1:length(rh)
                fprintf(fid,'%-12.8f\t',h1(i,j));
            end
            fprintf(fid,'\n');
        end
        fclose(fid);

        filename = 'h2data.m';
        fid=fopen(filename,'w');
        for i=1:lens
            for j=1:length(rh)
                fprintf(fid,'%12.8f\t',h2(i,j));
            end
            fprintf(fid,'\n');
        end
        fclose(fid);
        
        filename = 'y_data.m';
        fid=fopen(filename,'w');
        for i=1:length(rh)
            fprintf(fid,'%12.8f\n',rh(i));
        end
        fclose(fid);
        
    case 2 % headtype = 2 represents real head
        disp('headtype = 2: real head')

        y_coord = head_prop(11) + union(-linspace(0, head_prop(14), num_steps),linspace(0, head_prop(14), num_steps)); % so that the head can either be in front or behind island
        x_coord = head_prop(12) + head_prop(8); % realheadposition + head offset
        
        hx_av = zeros(length(y_coord),length(x_coord ),length(sv));
        hy_av = zeros(length(y_coord),length(x_coord ),length(sv));
        hz_av = zeros(length(y_coord),length(x_coord ),length(sv));

        head_prop(2) = 1; % to get field values in line with original data
        muo = 4*pi*1e-7; %in SI units

        lenx = length(x_coord);
        leny = length(y_coord);
        lens = length(sv);
        
        for k=1:lens
            islandgeo_prop(2) = a(k);
            islandgeo_prop(3) = b(k);
            islandgeo_prop(4) = t(k);
            islandgeo_prop(5) = alpha(k);
            
            for i=1:lenx
                for j =1:leny
                    % y coordinate is down track, x is cross track
                    head_prop(7) = head_prop(11) + islandgeo_prop(6) - y_coord(j); % since y_coord = head_prop(11) + islandgeo_prop(6) - head_prop(7) in real_headfield_av_vec
                    head_prop(8) = head_prop(12) + islandgeo_prop(7) - x_coord(i); % since x_coord = head_prop(12) + islandgeo_prop(7) - head_prop(8) in real_headfield_av_vec

                    % field converted to Tesla to speed up computation of volume average
                    % [hx_av(j,i) hy_av(j,i) hz_av(j,i)] = real_headfield_av_vec(realhead_pos_prop, realhead_field_prop*muo, head_prop, islandgeo_prop, tol_prop, step_vect)
                    [hx_av(j, i, k) hy_av(j, i, k) hz_av(j, i, k)] = real_headfield_av_vec(realhead_pos_prop, realhead_field_prop*muo, head_prop, islandgeo_prop, tol_prop, step_vect)
                end
            end
        end

        % converting to A/m
        hx_av = hx_av/muo;
        hy_av = hy_av/muo;
        hz_av = hz_av/muo;

        % will need to reshape the data after loading to memory
        
        filename = 'hx_av_data.m';
        fid=fopen(filename,'w');
        
        for k=1:lens
            for i=1:lenx
                for j=1:leny
                    fprintf(fid,'%-12.8f\t',hx_av(j, i, k)); % since we are printing a line horizontally
                end
                fprintf(fid,'\n');
            end
        end
        fclose(fid);

        filename = 'hy_av_data.m';
        fid=fopen(filename,'w');
        
        for k=1:lens
            for i=1:lenx
                for j=1:leny
                    fprintf(fid,'%-12.8f\t',hy_av(j, i, k)); % since we are printing a line horizontally
                end
                fprintf(fid,'\n');
            end
        end
        fclose(fid);

        filename = 'hz_av_data.m';
        fid=fopen(filename,'w');
        
        for k=1:lens
            for i=1:lenx
                for j=1:leny
                    fprintf(fid,'%-12.8f\t',hz_av(j, i, k)); % since we are printing a line horizontally
                end
                fprintf(fid,'\n');
            end
        end
        fclose(fid);
          
        filename = 'x_data.m';
        fid=fopen(filename,'w');
        for i=1:lenx
            fprintf(fid,'%12.8f\n',x_coord(i));
        end
        fclose(fid);

        filename = 'y_data.m';
        fid=fopen(filename,'w');
        for i=1:leny
            fprintf(fid,'%12.8f\n',y_coord(i));
        end
        fclose(fid);
        
    otherwise
        disp('unknown headtype')
end

filename = 's_data.m';
fid=fopen(filename,'w');
for i=1:lens
    fprintf(fid,'%12.8f\n',sv(i));
end
fclose(fid);
 end
