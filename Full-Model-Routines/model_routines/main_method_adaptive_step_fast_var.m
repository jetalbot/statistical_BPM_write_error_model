%% SWITCHING PROBABILITY SEARCH
% Finds the switching probability and the ideal next time step size to find
% the next value for the switching probability. Ideally this skips some areas
% and focuses on areas of interest.

% This function only returns switching probabilities for cases when there
% is a variation between islands.
%%

function prob_switch_data = main_method_adaptive_step_fast_var(var_prop, var_tol, tol_prop, tptiny_prop, adapt_prop, tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s)

    
    % Set initial values
    delta_xmin = adapt_prop(1); % smallest step size
    delta_xmax = adapt_prop(2); % largest step size
    delta_x = adapt_prop(3); %step size, first set to largest step size
    min_switch_prob = adapt_prop(4); % minimum switching probability
    switch_prob_low = adapt_prop(5); % low tolerance for switching probability
    switch_prob_high = adapt_prop(6); % high tolerance for switching probability

    rel_tol = var_tol(1); % relative tolerance for integration
    abs_tol = var_tol(2); % absolute tolerance for integration
    
    x_vec = []; %empty matrix to fill with the steps, x_vec = tsw
    prob_switch = [];%empty matrix to fill with the switching probability
    tsw= 0; % starting step size 
    exit_cnd = 1; % exit condition for step loops, if exit_cnd = 0 then loop stops
    i = 0; %loop count  
    n =0.2; % scales maximum possible step size in some regions
    
varparameter = var_prop(3);
s1 = var_prop(7);
s2 = var_prop(8);
w = var_prop(5);
dist_type = var_prop(6);

% Switching for variations, if no variations then use main_method_adaptive_step_fast_novar.m 
    switch varparameter
        case 1
            disp('shape variations results')
            filename = 'prob_switch_shape_var_temp.m';
            
        case 2
            disp('size variations results')
            filename = 'prob_switch_size_var_temp.m';
            
        case 3
            disp('position variations results')
            filename = 'prob_switch_position_var_temp.m';
            
        case 4
            disp('intrinsic anisotropy variations results')
            filename = 'prob_switch_k1_var_temp.m';
        case 5
            disp('island field variations results')
            filename = 'prob_switch_islandfield_var_temp.m';    
        otherwise
            disp('variations on')
    end


    fid=fopen(filename,'w');
    
    while exit_cnd %if exit_cnd = 0 loop exits
        
        disp('i = '); disp(i);
        
% Switch for kind of integration to us
if size(islandmag_prop)<6
        switch tol_prop(6)
            case 0
            tsw
% Switch for distribution type, at the moment, only Gaussian and truncated
% Gaussian

               switch dist_type
                   case 0     %if using gaussian variation case 
                    pswitch_int = quad(@pswitchIntegrand_fast, s1, s2, abs_tol, [], var_prop, tol_prop, tptiny_prop, tsw, tperiod,n, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s);
                    pswitch_val = pswitch_int
                   case 1 %if variations using the truncated gaussian case
                    pswitch_int = quad(@pswitchIntegrand_fast_nongauss, s1, s2, abs_tol, [], var_prop, tol_prop, tptiny_prop, tsw,tperiod, w, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s);
                    pswitch_val = pswitch_int
               end
            case 1
% Switch for distribution type     
                switch dist_type 
                    case 0
                    pswitch_int = quadgk(@(s)pswitchIntegrand_fast(s, var_prop, tol_prop, tptiny_prop, tsw, tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s), s1, s2, 'RelTol', rel_tol, 'AbsTol', abs_tol);
                    pswitch_val = pswitch_int
                    case 1
                    pswitch_int = quadgk(@(s)pswitchIntegrand_fast_nongauss(s, var_prop, tol_prop, tptiny_prop, tsw, tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s), s1, s2, 'RelTol', rel_tol, 'AbsTol', abs_tol);
                    pswitch_val = pswitch_int
                end        
            otherwise
                disp('unknown quadrature option')
        end
else
        switch tol_prop(6)
            case 0
            tsw
% Switch for distribution type, at the moment, only Gaussian and truncated
% Gaussian

               switch dist_type
                   case 0     %if using gaussian variation case 
                    pswitch_int = quad(@pswitchIntegrand_fast_ecc, s1, s2, abs_tol, [], var_prop, tol_prop, tptiny_prop, tsw, tperiod,n, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s);
                    pswitch_val = pswitch_int
                   case 1 %if variations using the truncated gaussian case
                    pswitch_int = quad(@pswitchIntegrand_fast_nongauss_ecc, s1, s2, abs_tol, [], var_prop, tol_prop, tptiny_prop, tsw,tperiod, w, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s);
                    pswitch_val = pswitch_int
               end
            case 1
% Switch for distribution type     
                switch dist_type 
                    case 0
                    pswitch_int = quadgk(@(s)pswitchIntegrand_fast_ecc(s, var_prop, tol_prop, tptiny_prop, tsw, tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s), s1, s2, 'RelTol', rel_tol, 'AbsTol', abs_tol);
                    pswitch_val = pswitch_int
                    case 1
                    pswitch_int = quadgk(@(s)pswitchIntegrand_fast_nongauss_ecc(s, var_prop, tol_prop, tptiny_prop, tsw, tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s), s1, s2, 'RelTol', rel_tol, 'AbsTol', abs_tol);
                    pswitch_val = pswitch_int
                end        
            otherwise
                disp('unknown quadrature option')
        end            
end    
        
% Print values to file        
        fprintf(fid,'%g %g\n',[tsw pswitch_val]);

% Find size for next step.        
% Here pswitch > 1 - p_low so we advance in position and step size is
% greater or equal to the maximum step size.
        if ((pswitch_int > 1 - switch_prob_low) && (delta_x >= delta_xmin) && (tsw >= 0)) 
            disp('(pswitch_int > 1 - switch_prob_low) && (delta_x >= delta_xmin) && (tsw >= 0)')
            i=i+1;
            x_vec(i) = tsw;
            prob_switch(i) = pswitch_val; 
            pswitch_int
            % Advance with the maximum possible step
            tsw = tsw + delta_xmax; 
            tsw
% Switching probability equal to more than at 1-switch_prob_high and equal
% to or less than 1-switch_prob_low, and the step size is greater than the
% minimum stepsize
        elseif ((pswitch_int >= 1 - switch_prob_high)&& (pswitch_int <= 1 - switch_prob_low) && (delta_x > delta_xmin) && (tsw >= 0)) 
            disp('(pswitch_int>= 1 - switch_prob_high)&& (pswitch_int/norm_factor <= 1 - switch_prob_low) && (delta_x > delta_xmin) && (tsw > 0)')
            delta_x = delta_x/2;
            % Here we keep going back until we the step_length is less than the mimimum allowed
            tsw = tsw - delta_x; 
            % Ensuring we dont get big negative numbers
            tsw = max(0,tsw); 
            tsw
% Switching probability is more than or equal to 1-switch_prob_high and
% less than or equal to 1-switch_prob_low, and the step size is less than
% the minimum size
        elseif ((pswitch_int >= 1 - switch_prob_high)&& (pswitch_int <= 1 - switch_prob_low) && (delta_x <= delta_xmin) &&(tsw >= 0)) 
            disp('(pswitch_int >= 1 - switch_prob_high)&& (pswitch_intr <= 1 - switch_prob_low) && (delta_x <= delta_xmin) &&(tsw >= 0)')
            i=i+1;
            x_vec(i) = tsw;
            prob_switch(i) = pswitch_val;
            pswitch_int
            % Advance with the minimum possible step
            tsw = tsw + delta_xmin; 
            tsw
% Switching probability greater than or equal to switch_prob_high and
% switching probability is less than 1-switch_prob_high
        elseif ((pswitch_int >= switch_prob_high)&& (pswitch_int < 1 - switch_prob_high) && (tsw >= 0)) 
            disp('(pswitch_int >= switch_prob_high)&& (pswitch_int < 1 - switch_prob_high)&& (tsw >= 0)')
            i=i+1;
            x_vec(i) = tsw;
            prob_switch(i) = pswitch_val; 
            pswitch_int
            % Advance with the maximum possible step
            tsw = tsw + n*delta_xmax; 
            tsw
 % Switching probability is greater than or equal to switch_prob_low and
 % less than switch_prob_high. And the step size is greater than the
 % minimum step size.
        elseif ((pswitch_int >= switch_prob_low)&& (pswitch_int < switch_prob_high) && (delta_x > delta_xmin) && (tsw >= 0)) 
            disp('(pswitch_int >= switch_prob_low)&& (pswitch_int < switch_prob_high) && (delta_x > delta_xmin) && (tsw >= 0)')
            delta_x = delta_x/2;
            % Here we keep going back until we the step_length is less than the mimimum allowed
            tsw = tsw - delta_x; 
            % Ensuring we dont get big negative numbers
            tsw = max(0,tsw); 
            tsw
 % Switching probability is greater than or equal to switch_prob_low and
 % less than switch_prob_high. And the step size is less than or equal to
 % the minimum step size.
        elseif ((pswitch_int >= switch_prob_low)&& (pswitch_int < switch_prob_high) && (delta_x <= delta_xmin) && (tsw >= 0)) 
            disp('(pswitch_int >= switch_prob_low)&& (pswitch_int < switch_prob_high) && (delta_x <= delta_xmin) && (tsw >= 0')
            i=i+1;
            x_vec(i) = tsw;
            prob_switch(i) = pswitch_val;
            pswitch_int
            % Advance with the minimum possible step
            tsw = tsw + delta_xmin; 
            tsw
% Switching probability is less than switch_prob_low and greater than
% minimum switching probability allowed.
        elseif ((pswitch_int < switch_prob_low)&& (pswitch_int > min_switch_prob) && (tsw >= 0)) 
            disp('(pswitch_int < switch_prob_low)&& (pswitch_int > min_switch_prob) && (tsw >= 0)')
            i=i+1;
            x_vec(i) = tsw;
            prob_switch(i) = pswitch_val; 
            pswitch_int
            % Advance with the maximum possible step
            tsw = tsw + n*delta_xmax; 
            tsw
        elseif tsw==0
            disp('tsw==0')
            i=i+1;
            x_vec(i) = tsw;
            prob_switch(i) = pswitch_val;
            pswitch_int
            % Advance with the maximum possible step
            tsw = tsw + delta_xmax; 
            tsw
% Here pswich > 1 - p_low so we advance in position. Switching probability
% is greater than 1-switch_prob_low
          elseif pswitch_int > 1 - switch_prob_low   
            disp('pswitch_int > 1 - switch_prob_low')
            i=i+1;
            x_vec(i) = tsw;
            prob_switch(i) = pswitch_val;
            pswitch_int
            % Advance with the maximum possible step
            tsw = tsw + delta_xmax; 
            tsw
            
        elseif pswitch_int < min_switch_prob
            disp('pswitch_int < min_switch_prob')
            exit_cnd = 0;
        else
            disp('unknown condition')
            exit_cnd = 0;
        end
    end
    % Close temporary file
    fclose(fid); 
    % Sort step size and switching probabilty vectors
    prob_switch_data = sortrows([x_vec'  prob_switch']);
