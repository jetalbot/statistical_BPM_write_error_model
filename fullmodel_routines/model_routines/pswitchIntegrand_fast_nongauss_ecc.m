% Used to integrate a truncated Gaussian distribution of the switching
% probability. Swap out pz to use another piecewise distribution function,
% see pz.m for details. 
function out = pswitchIntegrand_fast_nongauss_ecc(s,var_prop, tol_prop, tptiny_prop, tsw,tperiod,n, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s)

out = pz(s,var_prop,n).*pswitch_fast_ecc(s, var_prop, tol_prop, tptiny_prop, tsw, tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s);

end