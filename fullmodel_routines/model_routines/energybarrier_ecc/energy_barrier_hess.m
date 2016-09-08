%% ENERGY BARRIER CALCULATION FOR EXCHANGE COUPLED COMPOSITE ISLANDS
% Calculates the energy barrier of a double layer island using the method
% of calculation given in (Kalezhi et al,2012). Energy of coupled composite
% islands in an applied magnetic field are calculated here using a two-spin
% approximation that assumes that the hard and soft layers of the islands
% are uniformly magnetised and that the change in magnetisation angle
% occurs at the interface.
% For a single layer the energy barrier can be expressed as:

%  $E = K*V*sin^2(theta_h) - mu_0*M*V*cos(theta - theta_H)$

% For two layers the total energy assuming in-plane rotation of the field
% is expressed as:

% $E = K_h*V_h*sin^2(theta_h) + K_s*V_s*sin^2(theta_s)
% J*A*cos(theta_s - theta_h) - mu_0*M_h*V_h*cos(theta_h - theta_H)-
% mu_0*M_s*V_s*cos(theta_s - theta_H)-
% mu_0*[(M_s*V_s*M_s*V_h)/(4*pi*d^3)]*[2*S_2*cos(theta_s)*cos(theta_h) -
% S_1*cos(theta_s)*cos(theta_h)]$

% Where K is the total anisotropy of the material (including shape and
% crystalline).
% V is the layer volume
% t is the thickness of the material
% A is the cross-sectional area
% d is the spacing between the centres of soft and hard layers
% M is the saturation magnetisation
% theta_h and theta_s are the polar angles of the magnetisation in the hard
% and soft layers respectively
% theta_H is the applied field polar angle
% H is the magnitude of the applied field
% J is the interlayer exchange coupling constant
% The subscripts _s and _h denote the soft and hard layers
% S_1 and S_2 represent factors that depend on the relative position and
% shape of the layers and are calculated elsewhere. 

%% ENERGY BARRIER CALCULATION
% The energy barrier is calculated here for given values of the applied
% field in the hard and soft layers, hp_h_vec and hp_s_vec. Applied field
% values for a given head position can be calculated from scaledbarrier.m
% that will generate hp_h_vec and hp_s_vec from an imported field of
% head field values and associated field angles; theta_H_s_vec and
% theta_H_s_vec. Alternatively these values can be input directly as an
% array. 

% $h_dip = M_s.*V_s./(4*pi*d.^3.*h_k_h)$

% Where:
% M_s = saturation magnetisation for soft layer, similarly M_h is the
% saturation magnetisation for the hard layer
% V_s = volume of soft layer
% d is the spacing between the centres of the layers
 
% $h_k_h = 2*k_h_c/(muo*m_h) + M_h.*(nxx_h.*cos(phih_h).^2 + nyy_h.*sin(phih_h).^2 - nzz_h)$
% Where nxx_h nyy_h and nzz_h are demagnetisation factors that can be
% calculated in demagfactors.m from the island properties

% $k = (K_s.*V_s)/(K_h.*V_h)$
% $m = (M_s.*V_s).(M_h.*V_h)
%

% The ENERGY MINIMA are calculated as:
% $E' = (1/2)*(sin(theta_h)^2 + (K_s.*V_s)/(K_h.*V_h)*sin(theta_s)^2 - hp_h*cos(theta_h - theta_H_h) 
% - (M_s*V_s)/(M_h*V_h)*hp_s*cos(theta_s - theta_H_s) 
% - j_exch.*cos(theta_s - theta_h) -
% M_s.*V_s./(4*pi*d.^3.*h_k_h)*(2*S_2*cos(theta_s).*cos(theta_h) - S_1.*sin(theta_s).*sin(theta_h))$

% The SADDLE POINT is calculated as:
% $E'(x) = [sin(theta_h)*cos(theta_h) + hp_h*sin(theta_h - theta_H_h) -
% j_exch*sin(theta_s)- theta_h) +
% M_s.*V_s./(4*pi*d.^3.*h_k_h)*(2*S_2*cos(theta_s)*sin(theta_h)) +
% S_1*sin(theta_s).*cos(theta_h))

% E'(y) = (K_s*V_s)/(K_h.*V_h)*sin(theta_s)*cos(theta_s) +
% M_s.*V_s./(M_h.*V_h)*hp_s*sin(theta_s - theta_H_s) +
% j_exch*sin(theta_s-theta_h) +
% h_dip*(2*S_2*sin(theta_s)*cos(theta_h)+
% S_1*cos(theta_s).*sin(theta_h))]

function [ebarrier1 ebarrier2]= energy_barrier_hess(hp_h_vec, hp_s_vec, theta_H_h_vec, theta_H_s_vec, j_exch, h_dip, k, m, s1_factor, s2_factor)
options = optimset('Display','off','TolFun',1e-15, 'TolX',1e-15);


%Set array sizes
ebarrier1 = zeros(size(theta_H_h_vec));
ebarrier2 = zeros(size(theta_H_h_vec));


tol = 1e-6; % the tolerance used between two angles (e.g. if theta_1 - theta_2 <= tol then they are the same)
stop_count = 10; %if the energy barrier is stuck and the minimum cannot be found within the stop_count then it exits

for i=1:length(theta_H_h_vec)
    %% First minimum
    % Begin searching for the first minimum of the function between 0 and
    % 0.
    
    theta_H_h = theta_H_h_vec(i);
    theta_H_s = theta_H_s_vec(i);
    hp_h = hp_h_vec(i);
    hp_s = hp_s_vec(i);

    [th_m1,~,exitflag] =fminsearch(@(x)energy_min(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor),[0,0],options);

    n_count = 0;
    while (exitflag~=1 && n_count<=stop_count)
        disp('energy barrier calculation: stuck at 1st minimum')
        [th_m1,~,exitflag] =fminsearch(@(x)energy_min(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor),th_m1,options);
        n_count = n_count +1;
    end

    %% Second minimum
    % Begin searching for the first minimum between of the function between
    % -pi and -pi.
    [th_m2,~,exitflag] =fminsearch(@(x)energy_min(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor),[-pi,-pi],options);
    
    n_count = 0;
    while (exitflag~=1 && n_count<=stop_count)
        disp('energy barrier calculation: stuck at 2nd minimum')
        [th_m2,~,exitflag] =fminsearch(@(x)energy_min(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor),th_m2,options);
        n_count = n_count +1;
    end

    %% Direct root for finding saddlepoint
    % The saddle point of the function is found from the psoitions of the
    % first and second minimum
    
    % th_m1 is the first minimum
    % th_m2 is the second minimum
 
    % Checks that the two angles are the same within the given tolerance if
    % they are then the saddle point is searched for at half the value of
    % the second minimu. 
    if (abs(th_m2(1)-th_m1(1))<= tol && abs(th_m2(2)-th_m1(2))<= tol)
        th_s0 = th_m2/2;
    else 
        th_s0 = (th_m1 + th_m2)/2;
        th_m1_0 = th_m1;
        th_m2_0 = th_m2;
    end
    
    % Searching for the saddle point using the starting location found from
    % above
    [th_s,~,exitflag] = fsolve(@(x)direct_root(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor),th_s0,options);  % Call optimizer
    
    n_count = 0;
    set = [1 2 3 -3];
    while ((~ismember(exitflag, set)) && n_count<=stop_count)
        disp('energy barrier calculation: stuck at saddle point')
        th_s0 = th_s;
        [th_s,~,exitflag] = fsolve(@(x)direct_root(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor),th_s0,options);  % Call optimizer
        n_count = n_count +1;
    end
    
    % Find energy barrier
    x(1) = th_s(1); % first minimum of the magnetisation angle theta
    x(2) = th_s(2); % second minimum of the magnetisation angle theta
    x(3) = hp_h; % applied field in the hard layer
    x(4) = hp_s; % applied field in the soft layer
    
    
    % If the determinant of the Hessian, delta_s, is greater than zero then
    % the first energy barrier is zero and the second is calculated from
    % the energy barrer calculation. Otherwise both energy barriers are
    % calculated frpm the energy barrier calculation. 
    % Finally if the energy barrier is calculated as a negative value it is
    % assumed to be equal to zero. 
    delta_s = hessian_det(x, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor);

    if delta_s >= 0 || (abs(th_m2(1)-th_m1(1))<= tol && abs(th_m2(2)-th_m1(2))<= tol)
        ebarrier1(i) = 0;
        ebarrier2(i) = energy(hp_h, hp_s, theta_H_h, theta_H_s, th_s(1), th_s(2), j_exch, h_dip, k, m, s1_factor, s2_factor) - energy(hp_h, hp_s, theta_H_h, theta_H_s, th_m2(1), th_m2(2), j_exch, h_dip, k, m, s1_factor, s2_factor);
    else
        ebarrier1(i) = energy(hp_h, hp_s, theta_H_h, theta_H_s, th_s(1), th_s(2), j_exch, h_dip, k, m, s1_factor, s2_factor) - energy(hp_h, hp_s, theta_H_h, theta_H_s, th_m1(1), th_m1(2), j_exch, h_dip, k, m, s1_factor, s2_factor);
        ebarrier2(i) = energy(hp_h, hp_s, theta_H_h, theta_H_s, th_s(1), th_s(2), j_exch, h_dip, k, m, s1_factor, s2_factor) - energy(hp_h, hp_s, theta_H_h, theta_H_s, th_m2(1), th_m2(2), j_exch, h_dip, k, m, s1_factor, s2_factor);
    end
    
    if ebarrier1(i)< 0
        disp('energy barrier negative')
        ebarrier1(i)=0;
    end
    

end
end

%% Calculating first and second minimum
function out = energy_min(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor)
theta_h = x(1); % magnetisation field angle in the hard layer
theta_s = x(2); % magnetisation field angle in the soft layer
out = (1/2)*(sin(theta_h).^2 + k.*sin(theta_s).^2) - hp_h.*cos(theta_h - theta_H_h) - m.*hp_s.*cos(theta_s - theta_H_s) - j_exch.*cos(theta_s -theta_h) - h_dip.*(2*s2_factor*cos(theta_s).*cos(theta_h) - s1_factor.*sin(theta_s).*sin(theta_h));
end

%% To find saddlepoint
% Use fsolve of the function to find the saddle point of the energy barrier
% theta_h = x(1), field angle for the hard layer
% theta_s = x(2), field angle for the soft layer
function out = direct_root(x, hp_h, hp_s, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor)
out = [sin(x(1))*cos(x(1)) + hp_h*sin(x(1)-theta_H_h) - j_exch*sin(x(2)-x(1)) + h_dip*(2*s2_factor*cos(x(2))*sin(x(1)) + s1_factor*sin(x(2)).*cos(x(1)));
    k*sin(x(2))*cos(x(2)) + m*hp_s*sin(x(2)-theta_H_s) + j_exch*sin(x(2)-x(1)) + h_dip*(2*s2_factor*sin(x(2))*cos(x(1))+ s1_factor*cos(x(2)).*sin(x(1)))];
end

%% Hessian determinant
% Returns the determinant of the energy barrier hessian matrix
function [delta e_x e_y e_xx e_yy e_xy] = hessian_det(x, theta_H_h, theta_H_s, j_exch, h_dip, k, m, s1_factor, s2_factor)

% x(1) = th_s(1); % first minimum of the magnetisation angle theta
% x(2) = th_s(2); % second minimum of the magnetisation angle theta
% x(3) = hp_h; % applied field in the hard layer
% x(4) = hp_s; % applied field in the soft layer

e_x = sin(x(1))*cos(x(1)) + x(3)*sin(x(1)-theta_H_h) - j_exch*sin(x(2)-x(1)) + h_dip*(2*s2_factor*cos(x(2))*sin(x(1)) + s1_factor*sin(x(2)).*cos(x(1)));
e_y = k*sin(x(2))*cos(x(2)) + m*x(4)*sin(x(2)-theta_H_s)+ j_exch*sin(x(2)-x(1))+ h_dip*(2*s2_factor*sin(x(2))*cos(x(1)) + s1_factor*cos(x(2)).*sin(x(1)));

e_xx = cos(x(1))^2 - sin(x(1))^2 + x(3)*cos(x(1) - theta_H_h) + j_exch*cos(x(2)-x(1)) + h_dip*(2*s2_factor*cos(x(2))*cos(x(1)) - s1_factor*sin(x(2)).*sin(x(1)));
e_yy = k*(cos(x(2))^2 - sin(x(2))^2) + m*x(4)*cos(x(2)- theta_H_s) + j_exch*cos(x(2)-x(1)) + h_dip*(2*s2_factor*cos(x(2))*cos(x(1))- s1_factor*sin(x(2)).*sin(x(1)));
e_xy  = -j_exch*cos(x(2) - x(1)) + h_dip*(-2*s2_factor*sin(x(2))*sin(x(1)) + s1_factor*cos(x(2))*cos(x(1)));
delta = e_xx*e_yy - e_xy*e_xy;
end

%% Calculate energy
function out = energy(hp_h, hp_s, theta_H_h, theta_H_s, theta_h, theta_s, j_exch, h_dip, k, m, s1_factor, s2_factor)
out = (1/2)*(sin(theta_h).^2 + k.*sin(theta_s).^2) - hp_h.*cos(theta_h - theta_H_h) - m.*hp_s.*cos(theta_s - theta_H_s) - j_exch.*cos(theta_s -theta_h)- h_dip.*(2*s2_factor*cos(theta_s).*cos(theta_h) - s1_factor*sin(theta_s).*sin(theta_h));
end



% % REFERENCES
% [1] Kalezhi, J. & Miles, J.J., 2011. An Energy Barrier Model for Write Errors in Exchange-Spring Patterned Media. Magnetics, IEEE Transactions on, 47(10), pp.2540-2543.
% [2] Kalezhi, J., Belle, B.D. & Miles, J.J., 2010. Dependence of Write-Window on Write Error Rates in Bit Patterned Media. Magnetics, IEEE Transactions on, 46(10), pp.3752-3759.
% [3] Kalezhi, J., Miles, J.J. & Belle, B.D., 2009. Dependence of Switching
% Fields on Island Shape in Bit Patterned Media. Magnetics, IEEE
% Transactions on, 45(10), pp.3531-3534.