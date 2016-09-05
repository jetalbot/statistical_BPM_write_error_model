%% FIND SWITCHING PROBABILITY
% Gets the energy barrier for either an ECC or single layer islands and
% finds the switching probability.
% The switching probability is given by:
% $p_{switch} = 1 - p_1(t) = 1 - exp(-int_{t'-0}^t v(t') dt')$
% In the variation case the switching probability is given by:
% $p_{switch} = 1 - (int p(a)*p_1(a,t) da)/(int p(a) da)
% This is derived initially from the case when there is uniaxial
% anisotropy (two possible magnetisation orientations). The number of
% islands who magnetisation is in one orientation varies with time
% according to:
% $dn_1/dt + (v_12 + v_21)n_1 = v_21*n$
% Where n represents the number of islands the n_1 is the number of islands
% whose magnetisation is in orientation 1. v_12 and v_21 represent the
% transition rates from orientation 1 to 2 and 2 to 1, respectively. The
% transition rates are given by:
% $v_12 = f0 * exp(-E_{1,barrier}/(k_B*T))$ 
% $v_21 = f0 * exp(-E_{2,barrier}/(k_B*T))$
% Where f0 is the attempt frequency, k_B is the Boltzmann constant and T is
% the absolute temperature. E_{1,barrier} and E_{2,barrier} represent the
% energy barriers for magnetisation escaping from orientations 1 and 2
% respectively. 
% For a large number of islands the equation can be transformed to describe
% the probability of not switching by dividing by the total number of
% islands, n. Then we get:
% $dp_1/dt + (v_12 + v_21)p_1 = v_21$
% Where p_1 = n_1/n, which is the probability of remaining in the
% orientation 1. This can be integrated to obtain the first equation. For
% more details see references. 
%%
function probswitch= pswitch_fast_ecc(s, var_prop, tol_prop, tptiny_prop, tsw, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s)
% Set initial values into arrays
m_h = islandmag_prop(2).*ones(size(s));
h_k_h = islandmag_prop(3).*ones(size(s));
nxx_h = islandmag_prop(4).*ones(size(s));
nyy_h = islandmag_prop(5).*ones(size(s));
nzz_h = islandmag_prop(6).*ones(size(s));

islandgeo = islandgeo_prop(1);
a = islandgeo_prop(2).*ones(size(s));
b = islandgeo_prop(3).*ones(size(s));
t_h = islandgeo_prop(4).*ones(size(s));
alpha = islandgeo_prop(5).*ones(size(s));
islandposition_d =  islandgeo_prop(6).*ones(size(s));
islandposition_a =  islandgeo_prop(7).*ones(size(s));
v_h = islandgeo_prop(8).*ones(size(s));
probswitch = zeros(size(s));
lens = length(s);

% Find initial values
islandgeo = islandgeo_prop(1); %islandgeometry set to 1,2,3 for different island type

vel = head_prop(9); % head velocity in nm/ns
islandposition_down = islandgeo_prop(6); % downtrack island position in nm
head_pos = head_prop(7); % head position, in nm

dtrack_sep = islandposition_down - head_pos; % downtrack separation 
demag_tol = tol_prop(2);
tperiod = tptiny_prop(5);

tpmax = tptiny_prop(1); % twait = twait./tperiod;   % waiting time in normalized units (units of time period)
tupper = (tpmax + dtrack_sep./(vel.*tperiod) - tsw).*ones(size(s)); % Set array size for steps 

attfreq = thermal_prop(3); % attempt frequency*writeattempts 
rel_tol = tol_prop(1); % relative tolerance for integration
abs_tol = tol_prop(4); % absolute tolerance for integration

%% CALCULATE VARIATIONS
switch var_prop(3)
    
    case 0
        disp('no variations')
    case 1
        disp('shape variations')
        % For shape variations, volume and thickness are fixed, a*b =ao*bo,
        % thus $a^2*beta = ao^2*betao => a = ao*sqrt(betao/beta)$, letting
        % $beta = s*betao => a = ao*sqrt(1/s)=ao/sqrt(s)$
        % b = beta*a = betao*s*a = bo*s*a/ao = bo*s*(ao/sqrt(s))/ao =
        % bo*sqrt(s)
        
        a = islandgeo_prop(2)./sqrt(s); % islandgeo_prop(2) =  semi major axis, a
        b = islandgeo_prop(3).*sqrt(s); % islandgeo_prop(3) =  semi major axis, b
        
        % Get demagnetisation factors
        if (interp_prop(1)==1) % interp_demag =1; get values from interpolated data
            [nxx_h nyy_h nzz_h] = demagfactors_interp(s, demag_data); 
        else
            [nxx_h nyy_h nzz_h] = demagfactors(a, b, t_h, alpha, islandgeo, demag_tol);
        end

    case 2
        disp('size variations')
        % For size variations, thickness is fixed, s = a*b/(ao*bo), thus a*b =
        % s*ao*bo =sqrt(s)ao*sqrt(s)bo. thus a =sqrt(s)ao, b =sqrt(s)bo

        a = islandgeo_prop(2).*sqrt(s);
        b = islandgeo_prop(3).*sqrt(s);
        v_h = islandgeo_prop(8).*s;

        % Get demagnetisation values
        if (interp_prop(1)==1) % interp_demag =1; get values from interpolated data
            [nxx_h nyy_h nzz_h] = demagfactors_interp(s, demag_data); 
        else
            [nxx_h nyy_h nzz_h] = demagfactors(a, b, t_h, alpha, islandgeo, demag_tol);
        end 

    case 3
        disp('position variations')
        jitter_down_or_cross = var_prop(6); % jitter_down_or_cross % if = 0 then downtrack jitter, if =1 then crosstrack jitter
        % Switch depending on of cross-track or down-track jitter
        switch jitter_down_or_cross
            case 0
                islandposition_d = islandgeo_prop(9).*s; % downtrackperiod*s
                tupper = tpmax + s - tsw; % tpmax + rho/(vel*tperiod) - tsw;
            case 1
                islandposition_a = islandgeo_prop(10).*s; % crosstrackperiod*s
                tupper = tpmax + s - tsw; % tpmax + rho/(vel*tperiod) - tsw;
            otherwise
                disp('Supposed to be jitter, but it is neither down nor cross track jitter.')
        end

    case 4
        disp('k1 variations')
        h_k_h = islandmag_prop(3).*s;
    case 5
        disp('island field variations')
        h_k_h = islandmag_prop(3).*s;
    case 6
        disp('ms variations')
        m_h = islandmag_prop(2).*s;
        h_k_h = islandmag_prop(3)./s;
    otherwise
        disp('Unknown variation parameter.')
end


%% FIND ENERGY BARRIERS
% Calculate the energy barriers, hence the switching probability can be
% calculates, i.e. if energy barrier disappears then island has switched. 

% In the following, the drag tester writing method is followed, i.e. before
% switching head, hg =0.

for i=1:lens

    % Updating island properties
    islandgeo_prop(2) = a(i); % island semi major axis in nm
    islandgeo_prop(3) = b(i); % island semi minor axis in nm
    islandgeo_prop(4) = t_h(i); % hard layer thickness in nm
    islandgeo_prop(5) = alpha(i); % to complete island prop, set to 1
    islandgeo_prop(6) = islandposition_d(i); % downtrack island position in nm
    islandgeo_prop(7) = islandposition_a(i); % cross track island position in nm
    islandgeo_prop(8) = v_h(i); % volume of the hard layer in nm^3
    
    islandmag_prop(2) = m_h(i); % magnetisation saturation of the hard layer in A/m
    islandmag_prop(3) = h_k_h(i); % effect head field in the hard layer, in A/m 
    % demag terms
    islandmag_prop(4) = nxx_h(i); 
    islandmag_prop(5) = nyy_h(i);
    islandmag_prop(6) = nzz_h(i);

    % Set initial switch time
    switchtimes = 0;
    %Set step number
    j=1;
    
    % Generate vectors for substeps within larger steps (tiny steps!)
    % tptiny_prop(2) = 5*(twaitsteps -1) + 1
    % tptiny_prop(3) = upper limit of substeps
    % tptiny_prop(4) = number of substeps
    % Use these to find where the energy barrier vanishes
    tptiny = linspace(0, tupper(i), tptiny_prop(2));
    tptiny_sub = linspace(0, tptiny_prop(3), tptiny_prop(4));
    % Combine valued of the two arrays generated
    tptiny = union(tptiny_sub, tptiny);

    % Set some more arrays
    lentptiny = length(tptiny);
    tp_dat = zeros(size(tptiny));
    eb1 = zeros(size(tptiny));
    eb2 = zeros(size(tptiny));
    
    
    % Now find the energy barrier
    % While switchtimes is less than 1 (basically exit condition) and the
    % the number of sub steps is less the number of loops, j.
    while (switchtimes < 1) && (j <= lentptiny) % (switchtimes < 1) && (tptiny(j) <= tupper(i))
        % If island is single layer
        if size(islandmag_prop)<=6
            [ebarrier1 ebarrier2] = scaledbarrier_sl(tptiny(j), var_prop, tol_prop, tsw, tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, h_data, x_data, y_data, s_data, s(i), h_data_h, h_data_s); % considering first barrier only!       
            tp_dat(j) = tptiny(j);
            eb1(j) = ebarrier1;
            eb2(j) = ebarrier2;
        % Else island is ECC   
        else
            [ebarrier1 ebarrier2 headfielddata] = scaledbarrier_ecc(tptiny(j), var_prop, tol_prop, tsw, tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandmag_prop, islandgeo_prop, interp_prop, h_data, x_data, y_data, s_data, s(i), h_data_h, h_data_s); % considering first barrier only!       
            tp_dat(j) = tptiny(j);
            eb1(j) = ebarrier1;
            eb2(j) = ebarrier2;
        end    
        % If we get this then island switched! i.e first barrier vanished
        if ((~isreal(ebarrier1) || ebarrier1==0) && switchtimes < 1)  
            disp('barrier complex');
            switchtimes = 1;
        end
        % Add another numnber to the loop count
        j=j+1;
    end

%% CALCULATING THE SWITCHING PROBABILITY
    if(switchtimes < 1) % island did not switch

        tprob_dat = exp(-eb1)+ exp(-eb2);
        
        % Interpolate the sub step array, tp_dat, to find a value tprob_dat
        % within the switching probability array tp. Integrate this to find
        % the switching probability. 
        switch tol_prop(5)
            case 0
                switchrateintegral = quad(@switchrateIntegrand_interp, 0, tupper(i), abs_tol,[], tp_dat, tprob_dat);
            case 1
                switchrateintegral = quadgk(@(tp)switchrateIntegrand_interp(tp, tp_dat, tprob_dat), 0, tupper(i),'RelTol', rel_tol,'AbsTol', abs_tol);
            otherwise
                disp('unknown quadrature option')
        end
        %Scale with the rest of the equation. 
        expont = attfreq*tperiod.*switchrateintegral;
        probswitch(i) = 1 - exp(-expont);
        
    else % island switched
        probswitch(i) = 1;
    end
    
        filename = strcat('headfielddata.m'); 
    fid=fopen(filename,'a+');
    fprintf(fid,'%f %f %f %f %f %f %f\n', [headfielddata(1) headfielddata(2) headfielddata(3) headfielddata(4) headfielddata(5) headfielddata(6) probswitch(i)]);
    fclose(fid);
end
end


%% ARRAY DETAILS
% tptiny_prop(1)  = twait = twait./tperiod;   % waiting time in normalized units (units of time period)
% tptiny_prop(2)  = tptinysteps = 1*(twaitsteps -1) + 1;  %5*(twaitsteps -1) + 1; % to have 5*twaitsteps, used to check whether energy barrier vanishes
% tptiny_prop(3)  = t_upper_sub = 5;
% tptiny_prop(4)  = tptiny_sub_steps = 101;% number of tiny steps

% thermal_prop(1) = temp % temperature 
% thermal_prop(2) = kb % boltzmann constant
% thermal_prop(3) = attfreq*write_attempts % attfreq, attempt frequency = f0=1000*1e9, write attempts on target islands

% islandgeo_prop(1)  = islandgeo % set to 1,2,3 for different island types
% islandgeo_prop(2)  = a %semi major axis
% islandgeo_prop(3)  = b %semi minor axis
% islandgeo_prop(4)  = t_h %thickness of hard layer
% islandgeo_prop(5)  = alpha % ratio of top semi major axis to bottom semi major axis
% islandgeo_prop(6)  = islandposition_a %downtrack position
% islandgeo_prop(7)  = islandposition_d %crosstrack position
% islandgeo_prop(8)  = v_h % volume of hard layer
% islandgeo_prop(9)  = downtrackperiod 
% islandgeo_prop(10) = crosstrackperiod
% islandgeo_prop(11) = v_s %volume of soft layer
% islandgeo_prop(12) = area
% islandgeo_prop(13) = t_s %thickness of soft layer

% var_prop(1) = mean_var % normalised mean k1
% var_prop(2) = sigma_par % normalised standard deviation
% var_prop(3) = varparameter % 0-3 for different variations
% var_prop(4) = jitter_down_or_cross % if = 0 then downtrack jitter, if =1 then crosstrack jitter


% head_prop(1) = headtype %switches between using the interpolation for Karlqvist when =1 or a real head when =2.
% head_prop(2) = hg head gap field in A/m:  set to -1A/m, the error function will scale the field amplitude: hdatagen uses hg=1 to get head fields similar to original
% head_prop(3) = phih %not used for this field type
% head_prop(4) = gapsize % gap size in nm
% head_prop(5) = polesize % pole size in nm  
% head_prop(6) = flyheight % fly height in nm
% head_prop(7) = headposition_d % initial down track position of head in nm
% head_prop(8) = headposition_a % initial cross track position of head in nm
% head_prop(9) = vel % head velocity in nm/ns
% head_prop(10) = tau  %in ns: Schabes: head field rise time (constant)
% head_prop(11) = realheadposition_d % given initial down track position of head in nm
% head_prop(12) = realheadposition_a % given initial cross track position of head in nm
% head_prop(13) = interlayer % given interlayer spacing in nm
% head_prop(14) = downtrack_travel %head down track travel distance in nm: <=150 -a for real head, to cover two islands if head start to switch above island
% 

% If single layer
% islandmag_prop(1) = muo %
% islandmag_prop(2) = ms % saturation magnetisation
% islandmag_prop(3) = hk = 2*K1/mu0*ms
% islandmag_prop(4) = nxx % demag factor
% islandmag_prop(5) = nyy %
% islandmag_prop(6) = nzz % 

% If ECC
% islandmag_prop(1)  = muo %
% islandmag_prop(2)  = m_h % saturation magnetisation for hard layer
% islandmag_prop(3)  = h_k_h = 2*K1/mu0*ms, for hard layer
% islandmag_prop(4)  = nxx_h % demag factor, for hard layer
% islandmag_prop(5)  = nyy_h % 
% islandmag_prop(6)  = nzz_h %
% islandmag_prop(7)  = nxx_s % demag factor, for soft layer
% islandmag_prop(8)  = nyy_s %
% islandmag_prop(9)  = nzz_s %
% islandmag_prop(10) = h_k_s %  2*K1/mu0*ms, for soft layer
% islandmag_prop(11) = m_s % saturation magnetisation for soft layer
% islandmag_prop(12) = j_xc % exchange coupling
% islandmag_prop(13) = s1_factor  
% islandmag_prop(14) = s2_factor