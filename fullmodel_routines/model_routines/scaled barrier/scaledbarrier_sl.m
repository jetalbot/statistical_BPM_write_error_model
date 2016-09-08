%% SCALED BARRIER CALCULATION

% Used by the function pswitch_fast.m to scale the inputs for the energy
% barrier and produce the energy barrier values in appropriate units for
% further calculations.

% Can be considerably slimmed down by using global variables, possibly
% removed. 
function [eb1 eb2] = scaledbarrier_sl(tp, var_prop, tol_prop, tsw, tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, h_data, x_data, y_data, s_data, s, h_data_h, h_data_s)
%% FIND ENERGY BARRIERS FROM APPROPRIATE FUNCTION



temp = thermal_prop(1);
kb = thermal_prop(2);  %*1e-23 Boltzmann constant
muo = islandmag_prop(1);


ms = islandmag_prop(2);
hk = islandmag_prop(3);
vol = islandgeo_prop(8);

% Find the head field values and the corresponding values for angle
[hav phih theta_H] = head_pos(tp, tsw, tperiod, var_prop, head_prop, islandgeo_prop,realhead_pos_prop, realhead_field_prop, tol_prop,step_vect, interp_prop, y_data, x_data, h_data, h_data_h, h_data_s, s_data, s);

%Get effective magnetic field

nxx = islandmag_prop(4);
nyy = islandmag_prop(5);
nzz = islandmag_prop(6);
hkeff = hk + ms.*(nxx.*cos(phih).^2 + nyy.*sin(phih).^2 - nzz);

%Calculate scaled magnetic field for energy barrier function
h = hav./hkeff;

%Find energy barrier for field and head field angle

[ebarrier1, ebarrier2] = energybarrier(h, theta_H);

% Ebarrier =muo*ms*hkeff*vol*1e-4*ebarrier1/(kb*temp);
eb1 = muo*ms*vol*(1e-4)*hkeff.*ebarrier1/(kb*temp); % since vol in 10^-27 and kb in 10^-23
eb2 = muo*ms*vol*(1e-4)*hkeff.*ebarrier2/(kb*temp);


end

%%ARRAY DETAILS
%head_prop(1)  = headtype %1= karlqvist, 2=real head
%head_prop(2)  = hg % head gap field
%head_prop(3)  = phih % not used
%head_prop(4)  = gapsize % gap size
%head_prop(5)  = polesize % pole size
%head_prop(6)  = flyheight 
%head_prop(7)  = headposition_d % initial downtrack head position
%head_prop(8)  = headposition_a % initial crosstrack head position
%head_prop(9)  = vel %velocity
%head_prop(10) = tau % headfield rise time
%head_prop(11) = realheadposition_d % given initial head position downtrack
%head_prop(12) = realheadposition_a % given initial head posiotn crosstrack
%head_prop(13) = interlayer % given interlayer spacing
%head_prop(14) = downtracktravel % travel distance for head downtrack

%thermal_prop(1) = temp % temperature 
%thermal_prop(2) = kb % boltzmann constant
%thermal_prop(3) = attfreq*write_attempts % attfreq, attempt frequency = f0=1000*1e9, write attempts on target islands

%islandmag_prop(1) = muo %
%islandmag_prop(2) = ms % saturation magnetisation
%islandmag_prop(3) = hk = 2*K1/mu0*ms
%islandmag_prop(4) = nxx % demag factor
%islandmag_prop(5) = nyy %
%islandmag_prop(6) = nzz %




