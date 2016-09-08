% Used for integrating the a Gaussian distribution of the switching
% probability. Where px is a piecewise function of the Gaussian, this can
% be changed to just a Gaussian function here. I think I left it this way
% to see if something was working...
function out = pswitchIntegrand_fast(s, var_prop, tol_prop, tptiny_prop, tsw, tperiod,n, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s)


out = px(s,var_prop,n).*pswitch_fast(s, var_prop, tol_prop, tptiny_prop, tsw, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s);

end