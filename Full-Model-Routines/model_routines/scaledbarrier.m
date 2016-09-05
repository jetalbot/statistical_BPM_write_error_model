%% SCALED BARRIER CALCULATION

% Used by the function pswitch_fast.m to scale the inputs for the energy
% barrier and produce the energy barrier values in appropriate units for
% further calculations.

% Can be considerably slimmed down by using global variables, possibly
% removed. 
function [eb1 eb2] = scaledbarrier(tp, var_prop, tol_prop, tsw, tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, h_data, x_data, y_data, s_data, s, h_data_h, h_data_s)
%% FIND ENERGY BARRIERS FROM APPROPRIATE FUNCTION



temp = thermal_prop(1);
kb = thermal_prop(2);  %*1e-23 Boltzmann constant
muo = islandmag_prop(1);
ms = islandmag_prop(2);
hk = islandmag_prop(3);
vol = islandgeo_prop(8);

nxx = islandmag_prop(4);
nyy = islandmag_prop(5);
nzz = islandmag_prop(6);

hkeff = hk + ms.*(nxx.*cos(phih).^2 + nyy.*sin(phih).^2 - nzz);

h = hav./hkeff;


[hav_h hav_s phih_h phih_s thetah_h thetah_s] = head_pos_ecc(tp, tperiod, var_prop, head_prop, islandgeo_prop, y_data, x_data, h_data, h_data_h, h_data_s);

if length(islandmag_prop)>6 % here we have ECC structure
    nxx_h = islandmag_prop(4);
    nyy_h = islandmag_prop(5);
    nzz_h = islandmag_prop(6);
    ms_hard = islandmag_prop(2);
    vol_hard = islandgeo_prop(8);

    nxx_s = islandmag_prop(7);
    nyy_s = islandmag_prop(8);
    nzz_s = islandmag_prop(9);

    h_k_s = islandmag_prop(10); % in A/m
    ms_soft = islandmag_prop(11);

    h_k_h = islandmag_prop(3) + ms_hard.*(nxx_h.*cos(phih_h).^2 + nyy_h.*sin(phih_h).^2 - nzz_h); % in A/m
    k_h = muo*ms*h_k_h/2;
    k_s_c = muo*ms_soft*h_k_s/2;
    k_s = k_s_c + (1/2)*muo*ms_soft.^2.*(nxx_s.*cos(phih_s).^2 + nyy_s.*sin(phih_s).^2 - nzz_s);

    vol_soft = islandgeo_prop(11);
    k =k_s.*vol_soft./(k_h.*vol_hard);
    area = islandgeo_prop(12);
    j_xc = islandmag_prop(12);

    j_exch = area.*j_xc./(2*k_h.*vol_hard).*1e9; % since area/vol =1/nm = 1e9 /m

    m = ms_soft.*vol_soft./(ms.*vol_hard);
    t_h = islandgeo_prop(4);
    t_s = islandgeo_prop(13);

    d = (t_h + t_s)/2; % distance between hard and soft layer centres in nm
    h_dip = ms_soft.*vol_soft./(4*pi*d.^3.*h_k_h);
    
    s1_factor = islandmag_prop(13);
    s2_factor = islandmag_prop(14);

    hp_h = hav_h./h_k_h;
    hp_s = hav_s./h_k_h;

    % [ebarrier1, ebarrier2] = energy_barrier_ecc(h, thetah, j_exch(1), h_dip(1), k(1), m(1));
    %[ebarrier1, ebarrier2]= energy_barrier_ecc(hp_h, hp_s, theta_h_h, theta_h_s, j_exch(1), h_dip(1), k(1), m(1))
    [ebarrier1, ebarrier2]= energy_barrier_hess(hp_h, hp_s, thetah_h, thetah_s, j_exch(1), h_dip(1), k(1), m(1), s1_factor(1), s2_factor(1));
    eb1 = 2*k_h.*vol_hard*(1e-4).*ebarrier1/(kb*temp); % since vol_hard in 10^-27 and kb in 10^-23
    eb2 = 2*k_h.*vol_hard*(1e-4).*ebarrier2/(kb*temp);   
else % here we have single domain structure
    [ebarrier1, ebarrier2] = energybarrier(h, thetah);
    % Ebarrier =muo*ms*hkeff*vol*1e-4*ebarrier1/(kb*temp);
    eb1 = muo*ms*vol*(1e-4)*hkeff.*ebarrier1/(kb*temp); % since vol in 10^-27 and kb in 10^-23
    eb2 = muo*ms*vol*(1e-4)*hkeff.*ebarrier2/(kb*temp);
end
%hav_prop = [hp_h hp_s thetah_h thetah_s];
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


%if single layer
%islandmag_prop(1) = muo %
%islandmag_prop(2) = ms % saturation magnetisation
%islandmag_prop(3) = hk = 2*K1/mu0*ms
%islandmag_prop(4) = nxx % demag factor
%islandmag_prop(5) = nyy %
%islandmag_prop(6) = nzz %

%if ECC
%islandmag_prop(1)  = muo %
%islandmag_prop(2)  = m_h % saturation magnetisation for hard layer
%islandmag_prop(3)  = h_k_h = 2*K1/mu0*ms, for hard layer
%islandmag_prop(4)  = nxx_h % demag factor, for hard layer
%islandmag_prop(5)  = nyy_h % 
%islandmag_prop(6)  = nzz_h %
%islandmag_prop(7)  = nxx_s % demag factor, for soft layer
%islandmag_prop(8)  = nyy_s %
%islandmag_prop(9)  = nzz_s %
%islandmag_prop(10) = h_k_s %  2*K1/mu0*ms, for soft layer
%islandmag_prop(11) = m_s % saturation magnetisation for soft layer
%islandmag_prop(12) = j_xc % exchange coupling
%islandmag_prop(13) = s1_factor  
%islandmag_prop(14) = s2_factor



