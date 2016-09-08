%% HEAD FIELD POSITION AND ANGLES FOR ECC
% The function finds the head field and head field angle for a particular
% point of an island over a period of time given by 'tiny' time steps
% within a larger time period. Both of these time steps are used to
% determine the position of the head.

% The function assumes an ECC structure for the island for single layer
% islands use head_pos.m 


function [hav_h hav_s havx_h havx_s phi_H_h phi_H_s theta_H_h theta_H_s] = head_pos_ecc(tp, tsw, tperiod, var_prop, head_prop, islandgeo_prop, realhead_pos_prop, realhead_field_prop, tol_prop, interp_prop, y_data, x_data, h_data, h_data_h, h_data_s, s_data, s)
%% HEAD AND HEADFIELD PROPERTIES
% Uses the velocity of the head and large time steps (x) and the small time
% steps (tp) within those larger steps to determine current head position. 

    vel = head_prop(9); %velocity of the head
    headpos_d = vel*tperiod.*tp + vel*tperiod.*tsw  + head_prop(7); % assume head is behind island: only head position downtrack

    tau = head_prop(10);
    tau = tau./tperiod; % in normalized units
    hg = head_prop(2); % head field in the gap
    hgn = hg.*erf(2*tp./tau);  % since head starts switching at tp = 0
    % hgn = hg.*cuspulse(tp, tau, 1); % last parameter, 1, is normalized time period

    % Set array sizes, length from time period (tp)
    lentp = length(tp);
    havx = zeros(size(tp));
    havy = zeros(size(tp));
    havz = zeros(size(tp));

    havx_h = zeros(size(tp));
    havy_h = zeros(size(tp));
    havz_h = zeros(size(tp));

    hav_h = zeros(size(tp));
    theta_H_h = zeros(size(tp));
    phi_H_h = zeros(size(tp));

    havx_s = zeros(size(tp));
    havy_s = zeros(size(tp));
    havz_s = zeros(size(tp));

    hav_s = zeros(size(tp));
    theta_H_s = zeros(size(tp));
    phi_H_s = zeros(size(tp));

%% FIND HEAD FIELD VALUES AND ANGLES
    
% Find the head field values and associated field angles for each of the small time steps.  
    
for i=1:lentp
    head_prop(2) = hgn(i); % resetting the head gap field
    head_prop(7) = headpos_d(i); % setting new head position. head island separation calculated in karlqfield_av, real_headfield_av_vec
    % Find the volume averaged head field values in x y and z for hard and
    % soft layers. 
    [havx(i) havy(i) havz(i) havx_h(i) havy_h(i) havz_h(i) havx_s(i) havy_s(i) havz_s(i)]= real_headfield_av_vec_interp(var_prop, head_prop, islandgeo_prop, h_data, x_data, y_data, s_data, s, h_data_h, h_data_s);


    hav_h(i) = sqrt(havx_h(i).^2 + havy_h(i).^2 + havz_h(i).^2).*sign(hgn(i));
    hav_s(i) = sqrt(havx_s(i).^2 + havy_s(i).^2 + havz_s(i).^2).*sign(hgn(i));
        
    % Finding the head field angle for the energy barrier calculation    
    hratioz_h = havz_h(i)./hav_h(i);
    hratioz_s = havz_s(i)./hav_s(i);
    hratioyx_h = havy_h(i)./havx_h(i);
    hratioyx_s = havy_s(i)./havx_s(i);
        
    % Find the head field angle in the high anisotropy (hard) layer of the island
    %
    if isnan(hratioz_h)
        theta_H_h(i) = 0;
    else
        theta_test = acos(hratioz_h); % need only acute angles for energy barrier calculation
        if theta_test >= pi/2
            theta_H_h(i) = pi - theta_test;
        else
            theta_H_h(i) = theta_test;
        end
    end
    
    
    if isnan(hratioyx_h)
        phi_H_h(i) = 0;
    else
        phi_H_h(i) = atan(hratioyx_h);
    end
    
    % Now for the lower anisotropy (soft) layer
    if isnan(hratioz_s)
        theta_H_s(i) = 0;
    else
        theta_test = acos(hratioz_s); % need only acute angles for energy barrier calculation
        if theta_test >= pi/2
            theta_H_s(i) = pi - theta_test;
        else
            theta_H_s(i) = theta_test;
        end
    end
    
    
    if isnan(hratioyx_s)
        phi_H_s(i) = 0;
    else
        phi_H_s(i) = atan(hratioyx_s);
    end
     
end
