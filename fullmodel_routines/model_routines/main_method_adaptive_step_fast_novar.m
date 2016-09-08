%% SWITCHING PROBABILITY SEARCH
% Finds the switching probability and the ideal next time step size to find
% the next value for the switching probability. Ideally this skips some areas
% and focuses on areas of interest.

% This function only returns switching probabilities for cases when there
% is NO variation between islands.
%%

function prob_switch_data = main_method_adaptive_step_fast_novar(var_prop, tol_prop, tptiny_prop, adapt_prop, tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s)
    
    % Set initial values
    delta_xmin = adapt_prop(1); % smallest step size
    delta_xmax = adapt_prop(2); % largest step size
    delta_x = adapt_prop(3); %step size, first set to largest step size
    min_switch_prob = adapt_prop(4); % minimum switching probability
    switch_prob_low = adapt_prop(5); % low tolerance for switching probability
    switch_prob_high = adapt_prop(6); % high tolerance for switching probability

    x_vec = []; %empty matrix to fill with the steps, x_vec = tsw
    prob_switch = [];%empty matrix to fill with the switching probability
    tsw= 0; % starting step size 
    exit_cnd = 1; % exit condition for step loops, if exit_cnd = 0 then loop stops
    i = 0; %loop count  
    n =0.2; % scales maximum possible step size in some regions
    
    disp('no variations results')
    s = 1;

    %Set up temporary file
    filename = 'prob_switch_temp.m';
    fid=fopen(filename,'w');

while exit_cnd %if exit_cnd = 0 loop exits
%  
%      for tsw = 1.566:0.001:1.5665
        disp('i = '); disp(i);      
        if size(islandmag_prop)<7
            pswitchval = pswitch_fast(s, var_prop, tol_prop, tptiny_prop, tsw, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s)
        else
            pswitchval = pswitch_fast_ecc(s, var_prop, tol_prop, tptiny_prop, tsw, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s)
        end
%      end
        % Print values to file            
        fprintf(fid,'%g %4.20f\n',[tsw pswitchval]); % will need to sort data afterwards
        
% Find size for next step.        
% Here pswitch > 1 - p_low so we advance in position and step size is
% greater or equal to the maximum step size.         
        if ((pswitchval > 1 - switch_prob_low) && (delta_x >= delta_xmin) && (tsw >= 0)) 
            disp('(pswitchval > 1 - switch_prob_low) && (delta_x >= delta_xmin) && (tsw >= 0)')
            i=i+1;
            x_vec(i) = tsw;
            prob_switch(i) = pswitchval;
            pswitchval
            % Advance with the maximum possible step
            tsw = tsw + delta_xmax; 
            tsw
% Switching probability equal to more than at 1-switch_prob_high and equal
% to or less than 1-switch_prob_low, and the step size is greater than the
% minimum stepsize
        elseif ((pswitchval >= 1 - switch_prob_high)&& (pswitchval <= 1 - switch_prob_low) && (delta_x > delta_xmin) && (tsw >= 0)) 
            disp('(pswitchval >= 1 - switch_prob_high)&& (pswitchval <= 1 - switch_prob_low) && (delta_x > delta_xmin) && (tsw > 0)')
            delta_x = delta_x/2;
             % Here we keep going back until we the step_length is less
             % than the mimimum allowed
            tsw = tsw - delta_x; 
            % Ensuring we dont get big negative numbers
            tsw = max(0,tsw); 
% Switching probability is more than or equal to 1-switch_prob_high and
% less than or equal to 1-switch_prob_low, and the step size is less than
% the minimum size           
        elseif ((pswitchval >= 1 - switch_prob_high)&& (pswitchval <= 1 - switch_prob_low) && (delta_x <= delta_xmin) &&(tsw >= 0)) 
            disp('(pswitchval >= 1 - switch_prob_high)&& (pswitchval <= 1 - switch_prob_low) && (delta_x <= delta_xmin) &&(tsw >= 0)')
            i=i+1;
            x_vec(i) = tsw;
            prob_switch(i) = pswitchval;
            pswitchval
            % Advance with the minimum possible step
            tsw = tsw + delta_xmin; %HERE HERE HERE!!!!
            tsw
% Switching probability greater than or equal to switch_prob_high and
% switching probability is less than 1-switch_prob_high            
        elseif ((pswitchval >= switch_prob_high)&& (pswitchval < 1 - switch_prob_high) && (tsw >= 0)) 
            disp('(pswitchval >= switch_prob_high)&& (pswitchval < 1 - switch_prob_high)&& (tsw >= 0)')
            i=i+1;
            x_vec(i) = tsw;
            prob_switch(i) = pswitchval;
            pswitchval
            % Advance with the maximum possible step
            tsw = tsw + n*delta_xmax; 
            tsw
 % Switching probability is greater than or equal to switch_prob_low and
 % less than switch_prob_high. And the step size is greater than the
 % minimum step size.
        elseif ((pswitchval >= switch_prob_low)&& (pswitchval < switch_prob_high) && (delta_x > delta_xmin) && (tsw >= 0)) 
            disp('(pswitchval >= switch_prob_low)&& (pswitchval < switch_prob_high) && (delta_x > delta_xmin) && (tsw >= 0)')
            delta_x = delta_x/2;
             % Here we keep going back until we the step_length is less
             % than the mimimum allowed
            tsw = tsw - delta_x;
            % Ensuring we dont get big negative numbers
            tsw = max(0,tsw); 
            tsw

 % Switching probability is greater than or equal to switch_prob_low and
 % less than switch_prob_high. And the step size is less than or equal to
 % the minimum step size.
        elseif ((pswitchval >= switch_prob_low)&& (pswitchval < switch_prob_high) && (delta_x <= delta_xmin) && (tsw >= 0)) %switching probability is greater than or equal to 0.000001 and less than 0.01. And the step size is less than or equal to the minimum step size. 

            disp('(pswitchval >= switch_prob_low)&& (pswitchval < switch_prob_high) && (delta_x <= delta_xmin) && (tsw >= 0')
            i=i+1;
            x_vec(i) = tsw;
            prob_switch(i) = pswitchval;
            pswitchval
            tsw = tsw +delta_xmin; %AND HERE HERE HERE!!!!! %delta_xmin; % advance with the minimum possible step
            tsw
% Switching probability is less than switch_prob_low and greater than
% minimum switching probability allowed.
        elseif ((pswitchval < switch_prob_low)&& (pswitchval > min_switch_prob) && (tsw >= 0)) % switching probability is less than 0.000001 and greater than 1e-9. 
            disp('(pswitchval < switch_prob_low)&& (pswitchval > min_switch_prob) && (tsw >= 0)')
            i=i+1;
            x_vec(i) = tsw;
            prob_switch(i) = pswitchval;
            pswitchval
            % Advance with the maximum possible step
            tsw = tsw + n*delta_xmax; 
            tsw

        elseif tsw==0
            disp('tsw==0')
            i=i+1;
            x_vec(i) = tsw;
            prob_switch(i) = pswitchval;
            pswitchval
            % Advance with the maximum possible step
            tsw = tsw + delta_xmax; 
            tsw
% Here pswich > 1 - p_low so we advance in position. Switching probability
% is greater than 1-switch_prob_low
        elseif (pswitchval > 1 - switch_prob_low)
            disp('pswitchval > 1 - switch_prob_low')
            i=i+1;
            x_vec(i) = tsw;
            prob_switch(i) = pswitchval;
            pswitchval
            % Advance with the maximum possible step
            tsw = tsw + delta_xmax; 
            tsw
            
        elseif (pswitchval < min_switch_prob) % switching probability is less than 1e-9. 
            disp('pswitchval < min_switch_prob')
            exit_cnd = 0;
        else
            disp('unknown condition')
            % exit_cnd = 0;
        end
end
    fclose(fid);
    prob_switch_data = sortrows([x_vec' prob_switch']);
