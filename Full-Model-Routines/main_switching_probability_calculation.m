%% A STATISTICAL MODEL OF WRITE ERRORS IN  BIT PATTERNED MEDIA

% % Calculates the switching probability for an array of nanoislands using
% imported data files for volume averaged head fields. These can be 2D or
% 3D data sets. Island properies (both geometrical and magnetic)and lattice
% size and spacing must be specified. The head field properties for the
% imported head field data must also be specified. 
% Variations between individual islands can be set, temperature can also be
% changed from default 300K, along with the number of write attempts and 
% the attempt frequency. 
% The routines to find coercivity, energy barrier and switching field can
% be used apart from the rest of the code. Their locations are listed below.  
% 

% % LIST OF SECTIONS IN SWITCHING PROBABILITY CALCULATION 

% CHOOSE LATTICE TYPE -> select a either RECTANGULAR or HEXAGONAL 
% LATTICE types and change their period downtrack and cross track

% CHOOSE ISLAND GEOMETRY -> select island geometry, can be either
% a TRUNCATED ELLIPTIC CONE, an ELLIPTIC CYLINDER or a PRISM. Also specify
% DIMENSIONS and ISLAND POSITION.

% SPECIFY MAGNETIC PROPERTIES -> specify magnetic properties of
% islands. SATURATION MAGNETISATION, CRYSTALLINE ANISOTOPY and if ECC
% EXCHANGE COUPLING. Comment/uncomment as appropriate when changing between
% ECC and single domain. 

% HEAD FIELD PROPERTIES -> specify initial known constants for head; 
% POLESIZE, GAPSIZE, FLYHEIGHT and set INITIAL REAL HEAD POSITION. 

% SPECIFY THERMAL PROPERTIES  -> specify TEMPERATURE, number of WRITE
% ATTEMPTS and ATTEMPT FREQUENCY. 

% GENERATE TIME PARAMETERS -> specify the time parameters to be used
% to calculate the energy barrier over. 

% SET VARIATION PARAMETERS -> specify VARIATION in the islands, if
% varparameter=0 all islands are IDENTICAL, =1 for SHAPE VARIATION, =2 for 
% SIZE VARIATION, =3 for POSITION JITTER, =4 for CRYSTALLINE ANISOTROPY VARIATION 

% INTERPOLATION OF VARIATION PARAMETERS -> interpolation of variation 
% parameters, can specify NUMBER OF POINTS FOR INTERPOLATION.

% CALCULATE SWITCHING PROBABILITIES -> specify uniform or adaptive steps, 
% loads HEAD DATA and runs MAIN_METHOD_ADAPTIVE_STEP_LASTEST_FAST 
% (or MAIN_METHOD_FAST if uniform steps) to find PROBABILITY SWITCH DATA

%% PREAMBLE
clear
clear all

path(path, 'model_routines');

path(path, 'islands_field');


path(path, 'model_routines\field'); % nonuniform head field
path(path, 'model_routines\field\bpm_head'); % nonuniform head field
path(path, 'model_routines\field\karlqvistfield'); % nonuniform head field

path(path, 'model_routines\energybarrier_ecc'); 
path(path, 'model_routines\energybarrier'); 

path(path, 'model_routines\write_simulations_routines');
path(path, 'model_routines\scaled barrier');

path(path, 'variables');
path(path, 'variables\demagfactors');
path(path, 'variables\demagfactors\conedemag');
path(path, 'variables\demagfactors\cyldemag');
path(path, 'variables\demagfactors\prismdemag');
path(path, 'variables\exact_magnetostatic_inter_cyl_fast');

%% CHOOSE LATTICE TYPE

% RECTANGULAR: Here islands are arranged in a rectangular lattice set down
% track and off track periods comment/uncomment the following three lines
% to activate/deactivate this lattice type

latticetype = 1; % 'rectangular';
downtrackperiod = 18; %18; %8; %12.7;   % down track period in nm:
crosstrackperiod = 18; %18; %20; %12.7;  % cross track period in nm
downtrack_travel = 290; % 290;  % head down track travel distance in nm: <=150 -a for real head, to cover two islands if head start to switch above island

% HEXAGONAL: Islands are arranged in a hexagonal lattice set down track and
% cross track periods Comment/Uncomment the following three lines to
% activate/deactivate this lattice type

% latticetype = 2; %'hexagonal';
% downtrackperiod = 13.6;  % down track period in nm
% crosstrackperiod = 11.8;  % cross track period in nm: downtrackperiod*sqrt(3)/2
% downtrack_travel = 290;


%% CHOOSE ISLAND GEOMETRY

% TRUNCATED ELLIPTIC CONE:
% Specify island parameters
% Comment/uncomment the following lines to activate/deactivate this
% geometry

% islandgeo = 1; % 1 represents truncated elliptic cone
% a = 12.5/2; %6.7/2; % semi major axis in nm
% b = a; % semi minor axis in nm
% t = 6.25; % thickness in nm
% alpha = 0.66; % ratio of top to bottom semi major axis
% vol= pi*t.*a.*b.*((alpha + 1/2).^2 + 3/4)/3;

% ELLIPTIC CYLINDER: Recommended by INSIC
% Specify island parameters
% Comment/uncomment the following lines to activate/deactivate this
% geometry

islandgeo = 2;
a = 8/2; %6.5/2; % semi major axis in nm : Simon's specification
b = a; % semi minor axis in nm
t = 6; % thickness in nm
alpha = 1; % any number will do, needed to complete island_prop
vol= pi*a.*b.*t;

% PRISM
% Specify island parameters
% Comment/uncomment the following lines to activate/deactivate this
% geometry

% islandgeo = 3;
% a = (downtrackperiod/sqrt(2))/2;   % half length in nm: 50% fill factor
% b = a;            % half width in nm
% t = 10;        % height in nm
% alpha = 1; % any number will do, otherwise will not be used
% vol = 4*a.*b.*t;

% Island position values
islandposition_d = downtrackperiod; %25;  % down track position in nm, island set to same as downtrack period initially it is assumed the island is set forward of the head target
islandposition_a = 0;  % cross track position in nm

% ECC properties
t_h = 4; % hard layer thickness in nm , t_h + t_s = t % if working with single layer, this should be equal to t!
t_s = 2; % soft layer thickness in nm % if working with singlelayer, this should be zero!

area = pi*a*b;
v_h = pi*a.*b.*t_h; % volume in nm^3
v_s = pi*a.*b.*t_s;

% Array containing island properties
islandgeo_prop = [islandgeo a b t_h alpha islandposition_d islandposition_a v_h downtrackperiod crosstrackperiod v_s area t_s];


%% SPECIFY ISLAND MAGNETIC PROPERTIES

muo = 4*pi*1e-7; % permeability of free space in SI units

% SINGLE DOMAIN ISLAND PROPERTIES: comment out for ECC
 
% demag_tol = 1e-6; % tolerance used in some integrals of demagnetising factors
% ms = 700e3; % saturation magnetisation in A/m:  700emu/cc INSIC recommendation
% %hk = 8.048e5; %6.6682e5; % in A/m from maximum head field gradient
% k1 =8.36e5; % in J/m^3 muo*ms*hk/2;  %  crystalline anisotropy constant in J/m^3: since hk = 2*k1./(muo*ms) => k1 = muo*ms*hk/2
% hk = 2*k1/(muo*ms);
%  
% [nxx nyy nzz] = demagfactors(a, b, t, alpha, islandgeo, demag_tol) %calculate demag factors, see ref [3].
%  
% islandmag_prop = [muo ms hk nxx nyy nzz]; % for single layer properties, make sure to uncomment islandmag_prop below for ECC
%  

% ECC ISLAND PROPERTIES: comment out for single domain

demag_tol = 1e-4; % tolerance used in some integrals of demagnetising factors
m_h = 1000e3;  % hard layer saturation magnetisation in A/m: % 400 emu/cm^3
m_s = 1400e3; %1400e3; % soft layer saturation magnetisation in A/m: % 1000 emu/cm^3

[nxx_h  nyy_h nzz_h] = demagfactors(a, b, t_h, alpha, islandgeo, demag_tol)
[nxx_s  nyy_s nzz_s] = demagfactors(a, b, t_s, alpha, islandgeo, demag_tol)

k_h_c = 1.16E+06; %0.1*13*1e6; % 0.1*13.5*1e6;   %  hard layer crystalline anisotropy in J/m^3   % = 15*1e6 erg/cm^3  %1 erg/(cm^3)=0.1joules/(m^3)
k_s_c = 0.1*1*1e6; %0.1*1*1e6;    % soft layer crystalline anisotropy in J/m^3   % = 1*1e6 erg/cm^3

phi_h = pi/2; % though not needed for a cylinder

j_xc = 5*1e-3; % exchange coupling in J/m^2   % = 5erg/cm^2   % 1erg/(cm^2)=0.001joules/(m^2)

% Since hkeff = 2*k1eff/(muo*ms)= hk + ms.*(nxx.*cos(phih).^2 +
% nyy.*sin(phih).^2 - nzz); k1eff= k1 + 0.5*muo*ms^2.*(nxx.*cos(phih).^2 +
% nyy.*sin(phih).^2 - nzz);

k_h = k_h_c + (1/2)*muo*m_h.^2.*(nxx_h.*cos(phi_h).^2 + nyy_h.*sin(phi_h).^2 - nzz_h);
k_s = k_s_c + (1/2)*muo*m_s.^2.*(nxx_s.*cos(phi_h).^2 + nyy_s.*sin(phi_h).^2 - nzz_s);

h_k_h = 2*k_h_c/(muo*m_h) % in A/m, the shape term added in the scaled barrier function
h_k_s = 2*k_s_c/(muo*m_s) % in A/m, the shape term added in the scaled barrier function

k =k_s.*v_s./(k_h.*v_h);
j_exch = area.*j_xc./(2*h_k_h.*v_h).*1e9 % since area/v_h =1/nm = 1e9 /m
m = m_s.*v_s/(m_h.*v_h);

d = (t_h + t_s)/2; % distance between moments in nm
h_dip = m_s.*v_s./(4*pi*d.^3.*h_k_h);

z = d;
r = 1e-20; %eps;
t_1 = t_h;
t_2 = t_s;

a_1 = a;
a_2 = a;
[s1_factor, s2_factor, s3_factor, s4_factor] = S_factors_fast(z, r, a_1, a_2, t_1, t_2)

% For ECC type, uncomment islandmag_prop containing ECC properties
islandmag_prop = [muo m_h h_k_h nxx_h nyy_h nzz_h nxx_s nyy_s nzz_s h_k_s, m_s, j_xc, s1_factor, s2_factor];


%% HEAD FIELD PROPERTIES

% Parameters for Karlqvist field
polesize = 11; % pole size in nm 
gapsize = 25;  % gap size in nm
flyheight = 5; % fly height in nm
phih = pi/2; % not used for this field type

% Parameters for the real read
x_vect =0; y_vect = 0; z_vect = 0;

hx_data = 0; hy_data = 0; hz_data = 0; % used in real head

realheadposition_d = 300; %300; % given initial down track position of head in nm
realheadposition_a = 300; %300; % given initial cross track position of head in nm
interlayer = 1; % given interlayer spacing in nm

% Specify gap field, pole size, gap size, fly height, velocity ...
headtype = 2; %switches between using the interpolation for Karlqvist when =1 or a real head when =2. 
hg = -1; %-1; %-0.8; %-1; %-0.664; %-1; % head gap field in A/m:  set to -1A/m, the error function will scale the field amplitude: hdatagen uses hg=1 to get head fields similar to original
vel = 20; %5e-4; %5e-6; %20; %%  % % velocity in nm/ns

headposition_d = 0; % initial down track position of head in nm
headposition_a = 0; % initial cross track position of head in nm, this to be updated for adjacent track islands or jitter in the cross track direction
tau = 0.06; %in ns: Schabes: head field rise time (constant)

% Array containing head properties
realhead_pos_prop = [x_vect; y_vect; z_vect];

realhead_field_prop = cat(4, hx_data, hy_data, hz_data);
clear hx_data hy_data hz_data

head_prop = [headtype hg phih gapsize polesize flyheight headposition_d headposition_a vel tau realheadposition_d realheadposition_a interlayer downtrack_travel];


%% SPECIFY THERMAL PROPERTIES
% These are boltzmann constant, temperature, attempt frequency

temp = 100;  % in Kelvin: INSIC recommendation
kb = 1.3806503;  %*1e-23 Boltzmann constant

% Here attempt frequency is multiplied by write attempts, since the
% formular predicts this

write_attempts = 1; %10^4; % write attempts on target islands
attfreq = 1000e9; %100; %1e3;  %100*1e9; %24.3954; % attempt frequency, fo = 1000*1e9, in Hz. 

% For applied field = 0, Brown, PR, 130, 5, 1963 (also in PRB, 78,064430
% (2008)) predicts, in the high energy barrier approximation, fo =
% alpha*gamma*(muo*hk)*(k1*v/(pi*kb*T))^(1/2)/(1 + alpha^2) gamma =
% 1.760859770*1e11 rad^-1 s^-1 T^-1, choosing alpha =0.1 In this case, fo =
% alpha*gamma*muo*hk*sqrt(k1*vol*1e-27/(pi*kb*temp*1e-23))/(1 +alpha^2)
% which gives 6.193970e+010 Hz = 61.9397 GHz

% gamma = 1.760859770*1e11;
% alpha = 0.1; %damping parameter
% attfreq = alpha*gamma*muo*hk*sqrt(k1*vol*1e-27/(pi*kb*temp*1e-23))/(1 +alpha^2);
% attfreq = attfreq/1e9; %in GHz

thermal_prop = [temp kb attfreq*write_attempts]; % multiply attempt frequency by write attempts

%% GENERATE TIME PARAMETERS

tperiod = downtrackperiod/vel;    % time period in ns
twaitd = downtrack_travel/vel;     % for real head, chosen so that the largest separation does not exceed that permitted by the real head field data
twait = twaitd./tperiod;  % waiting time in normalized units (units of time period)

twaitsteps = 51; %51; %101; %51; %601; %31; % 101 to have 100 steps, waiting time steps
tptinysteps = 5*(twaitsteps -1) + 1; %5*(twaitsteps -1) + 1; % to have 5*twaitsteps, used to check whether energy barrier vanishes
t_upper_sub = 5; %5;
tptiny_sub_steps = 101;% number of tiny steps
tptiny_prop = [twait tptinysteps t_upper_sub tptiny_sub_steps tperiod];

%% SET VARIATION PARAMETERS

% Default for jiiter
jitter_down_or_cross = 0; % by default position variations are downtrack

% Choose parameters to vary pswitch and demagdatagen and fielddatagen
% should be updated

w = 1.5; %for non-gaussian function, n*sigma = cutoff from gaussian, can be any number when using gaussian
dist_type = 0; %for gaussian variations use 0, for non-gaussian function (pz) use 1

% NO VARIATIONS
% Here islands are all identical

varparameter = 0;
sigma_par = 0; % normalised standard deviation
mean_par = 1;
s1 = 1; s2 = 1; % upper and lower limit of integration same

% SHAPE
% Here the island volume is fixed as the island in-plane ellipticity 
% (b/a) is varied. For other ways of varying shape, relevant sections in

% varparameter = 1;
% sigma_par = 0.075; % normalised standard deviation, i.e sigma_shape/mean shape
% mean_par = 1; % normalised mean shape (ellipticity)
% n1 = 10;    % number determining upper limit of integration, see x1
% n2 = 10;    % number determining lower limit of integration, see x2
% s1 = mean_par - n1*sigma_par;  % upper limit of intergration
% s2 = mean_par + n2*sigma_par;  % lower limit of intergration

% SIZE
% Here the island in-plane ellipticity (b/a=1) and thickness is fixed as
% the volume is varied. For other ways of varying size, relevant sections in
% pswitch and demagdatagen and fielddatagen should be updated

% varparameter = 2;
% sigma_par = 0.075; % normalised standard deviation, i.e sigma_size/mean size
% mean_par = 1; % normalised mean size (volume)
% n1 = 10;    % number determining upper limit of integration, see x1
% n2 = 10;    % number determining lower limit of integration, see x2
% s1 = mean_par - n1*sigma_par;  % upper limit of intergration
% s2 = mean_par + n2*sigma_par;  % lower limit of intergration

% POSITION:JITTER
% Here the island position is varied from the expected position

% varparameter = 3;
% sigma_par = 0.075; % normalised standard deviation, i.e sigma_position/downtrackperiod
% mean_par = 1; %0 % normalised mean position
% n1 = 10;    % number determining upper limit of integration, see x1
% n2 = 10;    % number determining lower limit of integration, see x2
% s1 = mean_par - n1*sigma_par;  % upper limit of intergration
% s2 = mean_par + n2*sigma_par;  % lower limit of intergration


%LEAVE COMMENTED UNLESS VARYING JITTER TO CROSSTRACK 
% jitter_down_or_cross = 1; %position variatians are crosstrack


% % CRYSTALLINE ANISOTROPY
% % Here the crystalline anisotropy constant, k1, is varied from the expected
% % value

% varparameter = 4;
% %[mean_para, sigma_par]= dipolefield (a, t, m_h)
% sigma_par = 0.01; %0.085; % normalised standard deviation, i.e sigma_k1/mean k1
% mean_par = 1; %1 % normalised mean k1
% n1 = 20;    % number determining upper limit of integration, see x1
% n2 = 20;    % number determining lower limit of integration, see x2
% s1 = mean_par - n1*sigma_par;  % upper limit of intergration
% s2 = mean_par + n2*sigma_par;  % lower limit of intergration

% ISLAND PATTERN
% The pattern of surrounding islands is varied, all variations calculated
% in island_field folder

% varparameter = 5;
% [mean_para, sigma_par]= dipolefield (a, t, m_h)
% % sigma_par = 0.075;
% mean_par = 0; 
% n1 = 10;    % number determining upper limit of integration, see x1
% n2 = 10;    % number determining lower limit of integration, see x2
% s1 = mean_par - n1*sigma_par;  % upper limit of intergration
% s2 = mean_par + n2*sigma_par;  % lower limit of intergration


var_prop = [mean_par sigma_par varparameter jitter_down_or_cross w dist_type s1 s2];


%% TOLERANCES FOR INTEGRATION

var_tol = [1e-8 1e-10]; % [1e-4 0]; %[1e-3 1e-8]; % [1e-4 0]; % [1e-6 1e-10];  % relative and abolute tolerance for integrating variation parameters: default [1e-6 1e-10] to disable set to zero
pswitch_rel_tol = 1e-6; %1e-4;  % switching probability relative tolerance: default 1e-6
pswitch_abs_tol = 1e-6; %0; %1e-12; %0; % 1e-10; % switching probability absolute tolerance: default 1e-10: to disable set to 0
headfield_tol = 1e-6;    % used in averaging head field over island volume
step_vect = [0 0; 0 0; 0 0]; % [6 6; 6 6; 6 6]; % range of subset of the head field matrix for interpolation, all zeros means use all data
pswitch_quadrature_option = 1; % choose 0 for quad, 1 for quadgk
var_quadrature_option = 0; % choose 0 for quad, 1 for quadgk
tol_prop = [pswitch_rel_tol demag_tol headfield_tol pswitch_abs_tol pswitch_quadrature_option var_quadrature_option];


%% INTERPOLATION OF VARIATION PARAMETERS

% Initialising
interp_demag = 0; % for no interpolation in demagnetisation factors
interp_field = 0; % for no interpolation of volume averaged head fields
demag_data = 0; % for no interpolation of demagnetisation factors
h_data = 0; % for no interpolation of volume averaged head fields
x_data = 1; % for no interpolation of volume averaged head fields
y_data = 1; % for no interpolation of volume averaged head fields
s_data = 1; % for no interpolation of volume averaged head fields

% Uncomment the following lines to interpolate demag and field values. This speeds up
% computations remarkably!

% % demag field interpolation
% % Note:interpolation is NOT necessary for the prism case since expression is analytic
% 
npts =50; % number of points for interpolation
demagdatagen(a, b, t_h, alpha, islandgeo, varparameter, s1, s2, npts, demag_tol); % considering the hard layer only
demag_data = load('demagdata.m'); 
interp_demag = 1;



%% OBTAIN HEAD FIELD DATA INTERNAL
% Generating field values, i.e. h1, h2 for each size and various field
% locations for interpolation, uncommment if generating head field data
% within this .m file, if obtaining from another data file comment out. 

% snpts = npts; %50;
% num_steps =2;
% headfielddatagen(realhead_pos_prop, realhead_field_prop, head_prop, islandgeo_prop, tperiod, tptiny_prop, var_prop, tol_prop, step_vect, s1, s2, snpts, num_steps)

% % switch headtype
% %     case 1 % headtype =1, karlqvist type
% %         h_data = cat(3,load('h1data.m'), load('h2data.m'));
% %         y_data = load('y_data.m');
% %         s_data = load('s_data.m');
% %     case 2 % headtye = 2, real head
% %         x_data = load('x_data.m');
% %         y_data = load('y_data.m');
% %         s_data = load('s_data.m');
% % 
% %         hx_avdata = load('hx_av_data.m');
% %         hy_avdata = load('hy_av_data.m');
% %         hz_avdata = load('hz_av_data.m');
% % 
% %         if (length(s_data)==1)
% %             h_data = cat(3, hx_avdata, hy_avdata, hz_avdata);
% %         else % reshape head field data
% %             [hx_avdata hy_avdata hz_avdata] = reshape_headfield(x_data, y_data, s_data, hx_avdata, hy_avdata, hz_avdata);
% %             h_data = cat(4, hx_avdata, hy_avdata, hz_avdata);
% %         end
% %     otherwise
% %         disp('unknown head type')
% % end

interp_field = 1; 
interp_prop = [interp_demag interp_field];

%% CALCULATE SWITCHING PROBABILITIES

% Here we consider the head at various cross track positions



% OBTAINING HEAD DATA DATA EXTERNAL 
% Obtaining head data from file instead of from within this .m file

% hx_av_2d_data = load('.\karlqvist_single_pole_reflected_in _soft_underlayer\karlqvist_hx_data.m');
% hy_av_2d_data = load('.\karlqvist_single_pole_reflected_in _soft_underlayer\karlqvist_hy_data.m');
% hz_av_2d_data = load('.\karlqvist_single_pole_reflected_in _soft_underlayer\karlqvist_hz_data.m');
 
% %hard layer data sets
% hx_av_2d_h_data = load('.\karlqvist_single_pole_reflected_in _soft_underlayer\karlqvist_hx_data.m'); % set these to single layer ones if working with single layer islands
% hy_av_2d_h_data = load('.\karlqvist_single_pole_reflected_in _soft_underlayer\karlqvist_hy_data.m');
% hz_av_2d_h_data = load('.\karlqvist_single_pole_reflected_in _soft_underlayer\karlqvist_hz_data.m');
 
% %soft layer data sets
% hx_av_2d_s_data = load('.\karlqvist_single_pole_reflected_in _soft_underlayer\karlqvist_hx_data.m');
% hy_av_2d_s_data = load('.\karlqvist_single_pole_reflected_in _soft_underlayer\karlqvist_hy_data.m');
% hz_av_2d_s_data = load('.\karlqvist_single_pole_reflected_in _soft_underlayer\karlqvist_hz_data.m');
% 
% y_data = load('.\karlqvist_single_pole_reflected_in _soft_underlayer\karlqvist_y_data.m');

% hx_av_2d_data = load('.\bpmheadfield_11x22_x901y901\bpm11x22_hx_av_2d_data.m');
% hy_av_2d_data = load('.\bpmheadfield_11x22_x901y901\bpm11x22_hy_av_2d_data.m');
% hz_av_2d_data = load('.\bpmheadfield_11x22_x901y901\bpm11x22_hz_av_2d_data.m');
 
% hx_av_2d_h_data = load('.\bpmheadfield_11x22_x901y901\bpm11x22_hx_av_2d_h_data.m');
% hy_av_2d_h_data = load('.\bpmheadfield_11x22_x901y901\bpm11x22_hy_av_2d_h_data.m');
% hz_av_2d_h_data = load('.\bpmheadfield_11x22_x901y901\bpm11x22_hz_av_2d_h_data.m');
 
% hx_av_2d_s_data = load('.\bpmheadfield_11x22_x901y901\bpm11x22_hx_av_2d_s_data.m');
% hy_av_2d_s_data = load('.\bpmheadfield_11x22_x901y901\bpm11x22_hy_av_2d_s_data.m');
% hz_av_2d_s_data = load('.\bpmheadfield_11x22_x901y901\bpm11x22_hz_av_2d_s_data.m');
 
% y_data = load('.\bpmheadfield_11x22_x901y901\bpm11x22_y_data.m');

hx_av_2d_data = load('.\bpmheadfield_11x22(600x600)\bpm11x22_hx_av_2d_data.m');
hy_av_2d_data = load('.\bpmheadfield_11x22(600x600)\bpm11x22_hy_av_2d_data.m');
hz_av_2d_data = load('.\bpmheadfield_11x22(600x600)\bpm11x22_hz_av_2d_data.m');

hx_av_2d_h_data = load('.\bpmheadfield_11x22(600x600)\bpm11x22_hx_av_2d_h_data.m');
hy_av_2d_h_data = load('.\bpmheadfield_11x22(600x600)\bpm11x22_hy_av_2d_h_data.m');
hz_av_2d_h_data = load('.\bpmheadfield_11x22(600x600)\bpm11x22_hz_av_2d_h_data.m');

hx_av_2d_s_data = load('.\bpmheadfield_11x22(600x600)\bpm11x22_hx_av_2d_s_data.m');
hy_av_2d_s_data = load('.\bpmheadfield_11x22(600x600)\bpm11x22_hy_av_2d_s_data.m');
hz_av_2d_s_data = load('.\bpmheadfield_11x22(600x600)\bpm11x22_hz_av_2d_s_data.m');

y_data = load('.\bpmheadfield_11x22(600x600)\bpm11x22_y_data.m');

%% SPECIFY STEPS 

% UNIFORM STEP SIZE:
% Uncomment following lines for uniform step size

% step_cnd = 1;
% swtimesteps = 101; %301; % 61 to have 60 steps, number of head switching positions
% swtime =linspace(0, 10,swtimesteps); % head switch times in normalized units (units of time period)
% swtime = union(swtime,linspace(1, 4, 201)); % added more points between 2.5 and 4
% leswtime=length(swtime);
% swtime = 0; % head switches at t=0


% ADAPTIVE STEP SIZE
% Uncomment following lines for adaptive step size

step_cnd = 2;
delta_xmin =0.0001; %0.005; %0.004; % 0.005; % 0.0015; %0.0015; %0.015; % smallest step size
delta_xmax =0.01; %0.1; %0.04; % 0.01; 64*delta_xmin; % largest step size
delta_x = delta_xmax;
min_switch_prob = 1e-8; %1e-12;% minimum switching probability

switch_prob_low = 1e-8;
switch_prob_high = 1e-1;

adapt_prop = [delta_xmin delta_xmax delta_x min_switch_prob switch_prob_low switch_prob_high];

timestamp = now;

%% FIND SWITCHING PROBABILITY FROM GIVEN VALUES

% Create temporary file to save data as it is generated

filename = ('prob_switch_temp.m');
fid=fopen(filename,'w+'); % discard existing file if already exist
fclose(fid);


dir_name =strcat('off_track_prob_files_variation_parameter_',num2str(varparameter),'');
if exist(dir_name,'dir')==7
    disp('deleting files')
else
    mkdir(dir_name); %create the directory
end

head_pos = 0; %0:0.5:100; % head position in cross track direction, 0 implies that head is directly on main track, no offset
head_centre = 300; %300; % e.g. 0 for Karlvist field going from -40 to 40; to be set to the appropriate value , e.g 300 for a real head centred at x =300


for i=1:length(head_pos) 
    i

    [hx_avdata hy_avdata hz_avdata] = headfield_interp(y_data, y_data, hx_av_2d_data, hy_av_2d_data, hz_av_2d_data, head_centre, head_pos(i));
    [hx_av_s_data hy_av_s_data hz_av_s_data] =headfield_interp(y_data, y_data, hx_av_2d_s_data, hy_av_2d_s_data, hz_av_2d_s_data, head_centre, head_pos(i));
    [hx_av_h_data hy_av_h_data hz_av_h_data] =headfield_interp(y_data, y_data, hx_av_2d_h_data, hy_av_2d_h_data, hz_av_2d_h_data, head_centre, head_pos(i)); 
    
    h_data = cat(3, hx_avdata, hy_avdata, hz_avdata);
    h_data_h = cat(3, hx_av_h_data, hy_av_h_data, hz_av_h_data);
    h_data_s = cat(3, hx_av_s_data, hy_av_s_data, hz_av_s_data);  
    
    filename = strcat('h_data.m'); %outputs file with downtrack position in column 1 and switching probability in column 2
    fid=fopen(filename,'a+');
    fprintf(fid,'%f %f %f\n',h_data);
    fclose(fid);

%     for j=1:length(hx_avdata)
%     j    
%         if j>=1 && j<18
%             X = (j*10)
%             hx_island(j) = hx_avdata(j)+X;
%             hy_island(j) = hy_avdata(j)+X;
%             hz_island(j) = hz_avdata(j)+X;      
%         end
%         if j>18 && j<36
%             X = (j*10)+(j*10);
%             hx_island(j) = hx_avdata(j)+X;
%             hy_island(j) = hy_avdata(j)+X;
%             hz_island(j) = hz_avdata(j)+X;
%         end            
%         if j>36 && j<54
%             X = (j*10)+(j*10)+(j*10)
%             hx_island(j) = hx_avdata(j)+X;
%             hy_island(j) = hy_avdata(j)+X;
%             hz_island(j) = hz_avdata(j)+X;
%         end    
%         if j>54
%             X = (j*10)+(j*10)+(j*10)+(j*10)
%             hx_island(j) = hx_avdata(j)+X;
%             hy_island(j) = hy_avdata(j)+X;
%             hz_island(j) = hz_avdata(j)+X;
%         end
%         
%        
%     end

  
       

% 
% h_island_data = cat(3, hx_island, hy_island, hz_island);    
%      filename = strcat('h_island_data.m'); %outputs file with downtrack position in column 1 and switching probability in column 2
%      fid=fopen(filename,'a+');
%      fprintf(fid,'%f %f %f\n',h_island_data);
%      fclose(fid);   
%     
    
    if varparameter ==0
            disp('no variations')
            prob_switch_data1 = main_method_adaptive_step_fast_novar(var_prop, tol_prop, tptiny_prop, adapt_prop, tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s);
    elseif ismember (varparameter, [1 2 3 4 5])
            disp('variations')
            prob_switch_data1 = main_method_adaptive_step_fast_var(var_prop, var_tol, tol_prop, tptiny_prop, adapt_prop, tperiod, thermal_prop, realhead_pos_prop, realhead_field_prop, head_prop, step_vect, islandgeo_prop, islandmag_prop, interp_prop, demag_data, h_data, x_data, y_data, s_data, h_data_h, h_data_s);
    else
            disp('unknown variation parameter')
    end
    
    % SAVE SWITCHING PROBABILITIES   
    cd(dir_name);

    filename = strcat('prob_switch_data_sigma_',strcat(num2str(sigma_par*100),'_var_',num2str(dist_type),'_',num2str(timestamp),'.m')); %outputs file with downtrack position in column 1 and switching probability in column 2
    fid=fopen(filename,'a+');
    for j=1:length(prob_switch_data1(:,1))
        fprintf(fid,'%f %4.20f\n',[prob_switch_data1(j,1) prob_switch_data1(j,2)]);
        
    end
    
    fclose(fid);
    cd('..');  
end



%% REFERENCES
% [1] Kalezhi, J. & Miles, J.J., 2011. An Energy Barrier Model for Write Errors in Exchange-Spring Patterned Media. Magnetics, IEEE Transactions on, 47(10), pp.2540-2543.
% [2] Kalezhi, J., Belle, B.D. & Miles, J.J., 2010. Dependence of Write-Window on Write Error Rates in Bit Patterned Media. Magnetics, IEEE Transactions on, 46(10), pp.3752-3759.
% [3] Kalezhi, J., Miles, J.J. & Belle, B.D., 2009. Dependence of Switching Fields on Island Shape in Bit Patterned Media. Magnetics, IEEE Transactions on, 45(10), pp.3531-3534.


