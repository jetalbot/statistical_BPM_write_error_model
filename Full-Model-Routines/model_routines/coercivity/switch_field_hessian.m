function [h_switch_h h_switch_s] = switch_field_hessian(theta_H_h_vec, theta_H_s_vec, j_exch, h_dip, k, m, s1_factor, s2_factor)
options = optimset('Display','off','TolFun',1e-15, 'TolX',1e-15);
tol_h = 1e-8;
tol = 1e-8;

h_switch_h = zeros(size(theta_H_h_vec));
h_switch_s = zeros(size(theta_H_s_vec));
x = zeros(1,3);

stop_count = 10;

for i=1:length(theta_H_h_vec)
    exit_cnd = 1;
    theta_H_h = theta_H_h_vec(i)
    theta_H_s = theta_H_s_vec(i)

    hp_h = 0;
    hp_s = 0;
    delta_hp_h = 1e-4; %5e-4;
    delta_hp_s = 1e-4; %5e-4;

    th_m1_0 = [-1e-6,-1e-6];
    th_m2_0 = [-pi,-pi];

    while(exit_cnd)

        % finding first minimum
        [th_m1,fval,exitflag] =fminsearch(@(x)energy_min(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor), th_m1_0, options);

        n_count = 0;
        while (exitflag~=1 && n_count<=stop_count)
            disp('switching field calculation using hessian: stuck at 1st minimum')
            [th_m1,fval,exitflag] =fminsearch(@(x)energy_min(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor), th_m1, options);
            n_count = n_count + 1;
        end

        % finding second minimum
        [th_m2,fval,exitflag] =fminsearch(@(x)energy_min(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor), th_m2_0, options);

        n_count = 0;
        while (exitflag~=1 && n_count<=stop_count)
            disp('switching field calculation using hessian: stuck at 2nd minimum')
            [th_m2,fval,exitflag] =fminsearch(@(x)energy_min(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor),th_m2,options);
            n_count = n_count + 1;
        end

        % direct root for finding saddle point
        % [th_s,fval] = fsolve(@(x)direct_root(x, hp, theta_H, j_exch, h_dip, k, m),[-pi/2,-pi/2],options);  % Call optimizer

        if (abs(th_m2(1)-th_m1(1))<= tol && abs(th_m2(2)-th_m1(2))<= tol)
            th_s0 = th_m2/2;
        else
            th_s0 = (th_m1 + th_m2)/2;
            th_m1_0 = th_m1;
            th_m2_0 = th_m2;
        end

        [th_s,fval,exitflag] = fsolve(@(x)direct_root(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor),th_s0,options);  % Call optimizer

        n_count = 0;

        set = [1 2 3 -3];
        while ((~ismember(exitflag, set)) && n_count<=stop_count)
            disp('switching field calculation using hessian: stuck at saddle point')
            exitflag
            th_s
            [th_s,fval,exitflag] = fsolve(@(x)direct_root(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor),th_s,options);  % Call optimizer
            n_count = n_count + 1;
        end

        x(1) = th_s(1);
        x(2) = th_s(2);
        x(3) = hp_h;
        x(4) = hp_s;
        delta_s = hessian_det(x, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor);

        if(delta_hp_h<=tol_h) % && delta_s <=0)
            exit_cnd = 0;
        elseif(delta_s < 0)
            hp_h = hp_h - delta_hp_h;
            hp_s = hp_s - delta_hp_s;
            % theta_H_s = fieldangle_function(theta_H_h, theta_H_h_data, theta_H_s_data);
            % hp_s = field_function(hp_h, hp_h_data, hp_s_data);
        elseif(delta_s > 0 || (abs(th_m2(1)-th_m1(1))<= tol && abs(th_m2(2)-th_m1(2))<= tol))
            delta_hp_h =  0.5*delta_hp_h;
            hp_h = hp_h + delta_hp_h;
            delta_hp_s =  0.5*delta_hp_s;
            hp_s = hp_s + delta_hp_s;
            % theta_H_s = fieldangle_function(theta_H_h, theta_H_h_data, theta_H_s_data);
            % hp_s = field_function(hp_h, hp_h_data, hp_s_data);
        end

    end
    h_switch_h(i) = hp_h;
    h_switch_s(i) = hp_s;
end
end