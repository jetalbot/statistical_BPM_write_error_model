function probswitch = pswitch(s, var_prop, tol_prop, tptiny_prop, tsw, tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s)

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

dtrack_sep = islandgeo_prop(6)- head_prop(7); % only down track separation in this case
vel = head_prop(9);
% tptiny = linspace(0,tptiny_prop(1),tptiny_prop(2)); % linspace(0,twait,tptinysteps); % smaller tp intervals in normalized units (units of time period)
% tpmax = max(tptiny);
% lentptiny = length(tptiny);
tpmax = tptiny_prop(1);
tupper = (tpmax + dtrack_sep./(vel.*tperiod) - tsw).*ones(size(s));
demag_tol = tol_prop(2);

switch var_prop(3)
    case 0
        disp('no variations')
    case 1
        disp('shape variations')
        % for shape variations, volume and thickness fixed, a*b =ao*bo, thus
        % a^2*beta = ao^2*betao => a = ao*sqrt(betao/beta), letting beta = s*betao
        % => a = ao*sqrt(1/s)=ao/sqrt(s)
        % b = beta*a = betao*s*a = bo*s*a/ao = bo*s*(ao/sqrt(s))/ao = bo*sqrt(s)

        a = islandgeo_prop(2)./sqrt(s);
        b = islandgeo_prop(3).*sqrt(s);

        if (interp_prop(1)==1) % interp_demag =1;
            [nxx_h nyy_h nzz_h] = demagfactors_interp(s, demag_data); % get values from interpolated data
        else
            [nxx_h nyy_h nzz_h] = demagfactors(a, b, t_h, alpha, islandgeo, demag_tol);
        end

    case 2
        disp('size variations')
        % for size variations, thickness fixed, s = a*b/(ao*bo), thus a*b =
        % s*ao*bo =sqrt(s)ao*sqrt(s)bo. thus a =sqrt(s)ao, b =sqrt(s)bo

        a = islandgeo_prop(2).*sqrt(s);
        b = islandgeo_prop(3).*sqrt(s);
        v_h = islandgeo_prop(8).*s;

        if (interp_prop(1)==1) % interp_demag =1;
            [nxx_h nyy_h nzz_h] = demagfactors_interp(s, demag_data); % get values from interpolated data
        else
            [nxx_h nyy_h nzz_h] = demagfactors(a, b, t_h, alpha, islandgeo, demag_tol);
        end

    case 3
        % disp('position variations')
        switch var_prop(4)
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
        disp('ms variations')
        m_h = islandmag_prop(2).*s;
        h_k_h = islandmag_prop(3)./s;
    otherwise
        disp('Unknown variation parameter.')
end

% in the following, the drag tester writing method is followed, i.e. before
% switching head, hg =0.

probswitch = zeros(size(s));
lens = length(s);

attfreq = thermal_prop(3);
rel_tol = tol_prop(1);
abs_tol = tol_prop(4);

for i=1:lens
    % checking whether the energy barrier vanishes, i.e whether the island switches

    % updating island properties
    islandgeo_prop(2) = a(i);
    islandgeo_prop(3) = b(i);
    islandgeo_prop(4) = t_h(i);
    islandgeo_prop(5) = alpha(i);
    islandgeo_prop(6) = islandposition_d(i);
    islandgeo_prop(7) = islandposition_a(i);
    islandgeo_prop(8) = v_h(i);
    
    islandmag_prop(2) = m_h(i);
    islandmag_prop(3) = h_k_h(i);
    islandmag_prop(4) = nxx_h(i);
    islandmag_prop(5) = nyy_h(i);
    islandmag_prop(6) = nzz_h(i);

    switchtimes = 0;
    j=1;
    tptiny = linspace(0, tupper(i), tptiny_prop(2));
    lentptiny = length(tptiny);
    
    while (switchtimes < 1) && (j <= lentptiny) % (switchtimes < 1) && (tptiny(j) <= tupper(i))
        ebarrier = scaledbarrier(tptiny(j), var_prop, tol_prop, tsw, tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, h_data, x_data, y_data, s_data, s(i), h_data_h, h_data_s); % considering first barrier only!
        if ((~isreal(ebarrier) || ebarrier==0) && switchtimes < 1) % island switched i.e first barrier vanished
            % disp('barrier complex');
            switchtimes = 1;
        end
        j=j+1;
    end

    % calculating the switching probability
    if(switchtimes < 1) % island did not switch
        % switchrateintegral = quad(@switchrateIntegrand, 0, tupper(i), abs_tol,[], var_prop, tol_prop, tsw, tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, h_data, x_data, y_data, s_data, s(i)); 
       % tprob = switchrateIntegrand_interp(tp, tptiny_dat, tprob_dat) % verify how this works to improve speed
        switchrateintegral = quadgk(@(tp)switchrateIntegrand(tp,var_prop, tol_prop, tsw, tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, h_data, x_data, y_data, s_data, s(i), h_data_h, h_data_s), 0, tupper(i),'RelTol', rel_tol,'AbsTol', abs_tol);
        expont = attfreq*tperiod.*switchrateintegral;
        probswitch(i) = 1 - exp(-expont);
    else % island switched
        probswitch(i) = 1;
    end

end
end

