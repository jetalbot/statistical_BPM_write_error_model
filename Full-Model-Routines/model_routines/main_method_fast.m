function prob_switch_data = main_method_fast(varparameter, sigma_par, s1, s2, var_tol, var_prop, tol_prop, tptiny_prop, swtime, tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s)
%CALCULATES SWITCHING PROBABILITY IN UNIFORM STEPS

leswtime=length(swtime);
swtime  =reshape(swtime,leswtime,1); %so that we have a column

if varparameter==0
    disp('no variations results')
    prob_switch = zeros(size(swtime));
    s = 1;

    filename = 'prob_switch_temp.m';
    fid=fopen(filename,'a');
    for i=1:leswtime
        disp('i = '); disp(i);
        prob_switch(i) = pswitch_fast(s, var_prop, tol_prop, tptiny_prop, swtime(i), tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s);
        pswitchval = prob_switch(i)
        fprintf(fid,'%g %g\n',[swtime(i) prob_switch(i)]);
    end
    fclose(fid);
    prob_switch_data = [swtime prob_switch];

elseif varparameter==3
    disp('position variations results')
    norm_factor = sigma_par.*sqrt(2*pi)
    pswitchIntegral= zeros(size(swtime));

    rel_tol = var_tol(1);
    abs_tol = var_tol(2);
    
    filename = 'prob_switch_temp.m';
    fid=fopen(filename,'a');
    for i=1:leswtime
        disp('i = '); disp(i);
            
        switch tol_prop(6)
            case 0
                pswitchIntegral(i) = quad(@pswitchIntegrand_fast, s1, s2, abs_tol, [], var_prop, tol_prop, tptiny_prop, swtime(i), tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s);
                pswitch_val = pswitchIntegral(i)/norm_factor
                
            case 1
                pswitchIntegral(i) = quadgk(@(y)pswitchIntegrand_fast(y, var_prop, tol_prop, tptiny_prop, swtime(i), tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s), s1, s2, 'RelTol', rel_tol, 'AbsTol', abs_tol);
                pswitch_val = pswitchIntegral(i)/norm_factor
                
            otherwise
                disp('unknown quadrature option')
        end
        
        fprintf(fid,'%g %g\n',[swtime(i) pswitchIntegral(i)/norm_factor]);
    end
    fclose(fid);
    prob_switch = pswitchIntegral./norm_factor;
    prob_switch_data = [swtime prob_switch];

elseif ismember(varparameter, [1 2 4])
    disp('shape/size/k1 variations results')
    norm_factor = sigma_par.*(1 + erf(1./(sigma_par.*sqrt(2)))).*sqrt(2*pi)/2
    pswitchIntegral= zeros(size(swtime));

    rel_tol = var_tol(1);
    abs_tol = var_tol(2);

    filename = 'prob_switch_temp.m';
    fid=fopen(filename,'a');
    for i=1:leswtime
        disp('i = '); disp(i);
        
        switch tol_prop(6)
            case 0
                pswitchIntegral(i) = quad(@pswitchIntegrand_fast, s1, s2, abs_tol, [], var_prop, tol_prop, tptiny_prop, swtime(i), tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s);
                pswitch_val = pswitchIntegral(i)/norm_factor
                
            case 1
                pswitchIntegral(i) = quadgk(@(y)pswitchIntegrand_fast(y, var_prop, tol_prop, tptiny_prop, swtime(i), tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s), s1, s2, 'RelTol', rel_tol, 'AbsTol', abs_tol);
                pswitch_val = pswitchIntegral(i)/norm_factor
                
            otherwise
                disp('unknown quadrature option')
        end

        fprintf(fid,'%g %g\n',[swtime(i) pswitchIntegral(i)/norm_factor]);
    end
    fclose(fid);
    prob_switch = pswitchIntegral./norm_factor;
    prob_switch_data = [swtime prob_switch];
else
    disp('unknown variation parameter')
end

end