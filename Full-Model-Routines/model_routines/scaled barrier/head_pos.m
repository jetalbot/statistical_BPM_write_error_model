%% HEAD FIELD POSITION AND ANGLES FOR SINGLE LAYER
% This function finds the head field and head field angle values for a 'real head' at a
% particular point over an islands over a time period given by 'tiny' time
% steps with a larger time period. Both of these time steps are used to
% determine the position of the head and the angle of the field at that
% point. 

% This function assumes a single layer island for an island with an ECC
% structure used head_pos_ecc.m
%%
function [hav phi_H theta_H] = head_pos(tp, tsw, tperiod, var_prop, head_prop, islandgeo_prop,realhead_pos_prop, realhead_field_prop, tol_prop, step_vect, interp_prop, y_data, x_data, h_data, h_data_h, h_data_s, s_data, s)
 %% HEAD AND HEADFIELD PROPERTIES
 % Uses the velocity of the head and the large time steps (tsw) and the
 % small time steps within those larger steps (tp) to determine current
 % head position.
 
    vel = head_prop(9); % velocity of the head
    headpos_d = vel*tperiod.*tp + vel*tperiod.*tsw  + head_prop(7); % assume head is behind island: only head position downtrack

    tau = head_prop(10);
    tau = tau./tperiod; % in normalized units
    hg = head_prop(2); % head field in the gap
    hgn = hg.*erf(2*tp./tau);  % since head start switching at tp = 0
    % hgn = hg.*cuspulse(tp, tau, 1); % last parameter, 1, is normalized time period

    % Set array sizes, length from time period (tp)
    lentp = length(tp);
    havx = zeros(size(tp));
    havy = zeros(size(tp));
    havz = zeros(size(tp));

    hav = zeros(size(tp));
    theta_H = zeros(size(tp));
    phi_H = zeros(size(tp));

%% FIND HEAD FIELD VALUES AND ANGLES
%
% Now the head field values and associated field angles are found for each
% of the small time steps.
for i=1:lentp
    head_prop(2) = hgn(i); % reset the head gap field
    head_prop(7) = headpos_d(i); % setting new head position. head island separation calculated in karlqfield_av, real_headfield_av_vec
    
    %Now find the volume averaged head fields for x y and z
    [havx(i) havy(i) havz(i) havx_h(i) havy_h(i) havz_h(i) havx_s(i) havy_s(i) havz_s(i)]= real_headfield_av_vec_interp(var_prop, head_prop, islandgeo_prop, h_data, x_data, y_data, s_data, s, h_data_h, h_data_s);
    hav(i) = sqrt(havx(i).^2 + havy(i).^2 + havz(i).^2).*sign(hgn(i));

    %Find the field angles for the energy barrier calculation    
    hratio = havz(i)./hav(i);

        
    if isnan(hratio) % find the headfield angle for the island from the headfield values
        theta_H(i) = 0;
    else
        theta_test = acos(hratio); % need only acute angles for energy barrier calculation
        if theta_test >= pi/2
            theta_H(i) = pi - theta_test;
        else
            theta_H(i) = theta_test;
        end
    end
    
    hratio = havy(i)./havx(i);
    
    if isnan(hratio)
        phi_H(i) = 0;
    else
        phi_H(i) = atan(hratio);
    end
    
    
end
