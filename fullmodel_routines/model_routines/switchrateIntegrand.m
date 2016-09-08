function tprob = switchrateIntegrand(tp, var_prop, tol_prop, tsw, tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, h_data, x_data, y_data, s_data, s, h_data_h, h_data_s)

[eb1 eb2] = scaledbarrier(tp, var_prop, tol_prop, tsw, tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, h_data, x_data, y_data, s_data, s, h_data_h, h_data_s);
tprob = exp(-eb1)+ exp(-eb2);
end