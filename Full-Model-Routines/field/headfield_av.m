function [havx havy havz havx_h havy_h havz_h havx_s havy_s havz_s] = headfield_av(var_prop, realhead_pos_prop, realhead_field_prop, head_prop, islandgeo_prop, tol_prop, step_vect, interp_prop, h_data, x_data, y_data, s_data, s, h_data_h, h_data_s)

switch head_prop(1)
    case 1
%         disp('Karlqvist field')
        if (interp_prop(2)==1) % interp_field =1;
            [havx havy havz] = karlqfield_av_interp(var_prop, head_prop, islandgeo_prop, h_data, y_data, s_data, s);
        else
            [havx havy havz] = karlqfield_av(head_prop, islandgeo_prop, tol_prop);
        end

    case 2
%         disp('real head')
        if (interp_prop(2)==1) % interp_field =1;
            [havx havy havz havx_h havy_h havz_h havx_s havy_s havz_s]= real_headfield_av_vec_interp(var_prop, head_prop, islandgeo_prop, h_data, x_data, y_data, s_data, s, h_data_h, h_data_s);
        else
            % here head field not yet averaged over island volume
            [havx havy havz] = real_headfield_av_vec(realhead_pos_prop, realhead_field_prop, head_prop, islandgeo_prop, tol_prop, step_vect);
        end

    otherwise
        disp('Unknown variation parameter.')
end

end