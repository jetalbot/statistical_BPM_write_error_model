%% SET VARIABLES FOR ALL SWITCHING PROBABILITY CALCULATIONS
% Here all the geometric and magnetic properties of an island can be set
% along with the details of the head and the thermal properties. 
%%

% SET PATHS
path(path, 'demagfactors');
path(path, 'demagfactors\conedemag');
path(path, 'demagfactors\cyldemag');
path(path, 'demagfactors\prismdemag');

path(path, 'exact_magnetostatic_inter_cyl_fast');


%% CHOOSE LATTICE TYPE

% Islands can be set in either a RECTANGULAR or HEXAGONAL lattice. 

% Set lattice type, 1 = rectangular, 2 = hexagonal
latticetype = 1; % Here islands are arranged in a rectangular lattice. 

downtrackperiod = 18; %down-track period in nm
crosstrackperiod = 18; %cross-track period in nm

% Set the distance the head travels down track in nm. For a real head the
% distance should be <=150nm to cover two islands if the head starts to
% switch above the previously written island
downtrack_travel = 200; 

%% CHOOSE ISLAND GEOMETRY

% An island can be set to be a truncated elliptic cone, an elliptic
% cylinder or a prism.

% Specify island geometry, 1 = elliptic cone, 2 = elliptic cylinder, 3 =
% prism. 
islandgeo = 2; % The island is set to be an elliptic cylinder

% Set dimensions
a = 8/2; % Semi major axis in nm(for the prism this is the half length)
b = a; % Semi minor axis in nm(for the prism this is the half width)
t = 6; % Thickness or height in nm
alpha = 1; % Ratio of top to bottom semi major axes, used only for elliptic cone, otherwise set to 1

% Island position values. The island is set to same value as downtrack
% period initially it is assumed the island is set forward of the head
% target.
islandposition_d = downtrackperiod; % downtrack position in nm
islandposition_a = 0 % crosstrack position in nm, 0 means no offset

% ECC layer thicknesses. These must be set so that t = t_h + t_s
t_h = 4; % hard layer thickness in nm, if single layer set to the same value as t
t_s = 2; % soft layer thickness in nm, if single layer set to zero

area = pi*a*b; % island area
% Volume of each layer
v_h = pi*a.*b.*t_h 
v_s = pi*a.*b.*t_s

% Make array for all island properties

islandgeo_prop = [islandgeo a b t_h alpha islandposition_d islandposition_a v_h downtrackperiod  crosstrackperiod v_s area t_s];

%% SPECIFY ISLAND MAGNETIC PROPERTIES

muo = 4*pi*1e-7; % Permeability of free space in SI units
demag_tol = 1e-6 % Tolerance used for some of the integrals for the demagnetisation factors

% SINGLE DOMAIN ISLAND PROPERTIES
% Use these values for a single layer island, use the ECC properties below
% for multilayer. 

ms = 500*1000; % Saturation magnetisation in A/m, 700 emu/cc INSIC recommendation
k1 = 3.55e5; % Crystalline anisotropy in J/m^3,
hk = 2*k1/(muo*ms)

% Calculate the demagnetisation factors for a single layer island
[nxx nyy nzz] = demagfactors(a, b, t, alpha, islandgeo, demag_tol)

% ARRAY FOR SINGLE LAYER PROPERTIES, comment out to use ECC
% islandmag_prop = [muo ms hk nxx nyy nzz]

% ECC ISLAND PROPERTIES
% Use these values for a multilayer island, for single layer use the ones
% above. 

m_h = 700*1000; %Hard layer saturation magnetisation in A/m, 
m_s = 1400*1000; %Soft layer magnetisation saturation in A/m

% Calculate demagnetisation factors for hard and soft layers
[nxx_h nyy_h nzz_h] = demagfactors(a, b, t_h, alpha, islandgeo, demag_tol)
[nxx_s nyy_s nzz_s] = demagfactors(a, b, t_s, alpha, islandgeo, demag_tol)


k_h_c = 13*0.1*1e6; % Hard layer crystalline anisotropy in J/m^3
k_s_c = 1*0.1*1e6; % Soft layer crystalline anisotropy in J/m^3

phi_H = pi/2; % Angle used to normalise the anisotropy values, not needed for cylindrical island

j_xc = 5*1e-3; % Exchange coupling in J/m^2, usually set to approximately 5erg/cm^2

% Calculate the effective anistropies
% Since the effective headfield: 
% $H_k_eff = 2*K_1_eff/muo*ms = H_k + M_s*(N_xx*cos(phi_H)^2 + N_yy*cos(phi_H)^2 + N_zz)$ 
% Then: 
% $K_1_eff = K_1 + (1/2)*muo*M_s^2*(N_xx*cos(phi_H)^2 + N_yy*cos(phi_H)^2 + N_zz)$

k_h = k_h_c + (1/2)*muo*m_h.^2.*(nxx_h.*cos(phi_H).^2 + nyy_h.*sin(phi_H).^2 - nzz_h);
k_s = k_s_c + (1/2)*muo*m_s.^2.*(nxx_s.*cos(phi_H).^2 + nyy_s.*sin(phi_H).^2 - nzz_s);

% Calculate effective head field, the shape term is added in the scaled
% barrier function. 
h_k_h = 2*k_h_c/(muo*m_h); % Effective head field in A/m
h_k_s = 2*k_s_c/(muo*m_s); % Effective head field in A/m

% Scaled anisotropy to be used in energy barrier equations
k = k_s.*v_s./(k_h.*v_h);

% Scaled exchange coupling 
j_exch = area.*j_xc./(2*h_k_h.*v_h); % since area/v_h = 1/nm = 1e9 /m

% Scaled saturation magnetisation 
m = m_s.*v_s/(m_h.*v_h);

%Distance between moments in nm
d = (t_h + t_s)/2; 

h_dip = m_s.*v_s./(4*pi*d.^3.*h_k_h);


% Dipolar interaction, S factor values. S1 and S2 represent factors that
% depend on the relative position and shape of the layers.

%Coordinates for the calculations, such that position vector rho =
%(r*cos(theta_r), r*sin(theta_r), z), etc.
z = d; % 
r = 1e-20; % 
t_1 = t_h; %
t_2 = t_s;

a_1 = a;
a_2 = a;

[s1_factor s2_factor s3_factor s4_factor] = S_factors_fast(z, r, a_1, a_2, t_1, t_2);

% ARRAY FOR ECC PROPERTIES, comment out to use single layer
islandmag_prop = [muo m_h h_k_h nxx_h nyy_h nzz_h nxx_s nyy_s nzz_s h_k_s, m_s, j_xc, s1_factor, s2_factor];

%% HEAD FIELD PROPERTIES

% Parameters for field
polesize = 11; % pole size in nm
gapsize = 25; %gap size between pole and return pole in SUL in nm
flyheight = 5; % fly height of head in nm
phih = pi/2; % azimuthal angle
headtype = 2; % switches between Karlqvist = 1 and real head = 2
vel = 20; % head velocity in nm/ns 

% Head gap field, setting to -1A/m and allow the error function to scale
% the field amplitude: hdatagen uses hg = 1 to get head field similar to
% the original. The drag tester writing method is followed when finding the
% energy barrier, i.e. before switching head, hg =0.
hg = -1; % head gap field in A/m, 

% Parameters for the real head
x_vect = 0; y_vect = 0; z_vect = 0;

% Real head starting values
hx_data = 0; hy_data = 0; hz_data = 0;

realheadposition_d = 300; % given initial down track position of the head in nm
realheadposition_a = 300; % given initial cross track position of the head in nm
interlayer = 1; % given interlayer spacing in nm

headposition_d = 0; %initial down track position of head in nm
headposition_a = 0; %initial cross track position of head in nm
tau = 0.06; %head field rise time (constant), in ns

% Array containing head properties
realhead_pos_prop = [x_vect; y_vect; z_vect];

realhead_field_prop = cat(4, hx_data, hy_data, hz_data);
clear hx_data hy_data hz_data

head_prop = [headtype hg phih gapsize polesize flyheight headposition_d headposition_a vel tau realheadposition_d realheadposition_a interlayer downtrack_travel];

%% SPECIFY THERMAL PROPERTIES
% Set values for temperature, and attempt frequency

temp = 300;  % in Kelvin: INSIC recommendation
kb = 1.3806503*1e-23; %Boltzmann constant

% Here attempt frequency is multiplied by write attempts, since the
% formular predicts this

write_attempts = 1; % number of write attempts on target islands
attfreq = 100e9; % attempt frequency, fo = 1000*1e9, in Hz. 

thermal_prop = [temp kb attfreq*write_attempts]; % multiply attempt frequency by write attempts

%% GENERATE TIME PARAMETERS

tperiod = downtrackperiod/vel; % time period in ns
twaitd = downtrack_travel/vel; % for real head, chosen so that the largest separation does not exceed that permitted by the real head field data
twait = twaitd./tperiod;  % waiting time in normalized units (units of time period)

twaitsteps = 51; % 101 to have 100 steps, waiting time steps

% To check whether energy barrier vanishes: 5*(twaitsteps -1) + 1 to have
% 5*twaitsteps
tptinysteps = 5*(twaitsteps -1) + 1; t_upper_sub = 5; %5;
tptiny_sub_steps = 101;% number of tiny steps

% Array of time parameters
tptiny_prop = [twait tptinysteps t_upper_sub tptiny_sub_steps tperiod];

%% SET VARIATION PARAMETERS

% Default for jiiter
jitter_down_or_cross = 0; % by default position variations are downtrack

% CHOOSE DISTRIBUTION
% For gaussian variations, choose dist_type = 0 for truncated gaussian
% choose dist_type = 1. For more distributions update pz.m, or create a new
% piecewise function. 
dist_type = 0; % 

w = 1.5; % for non-gaussian function, n*sigma = cutoff from gaussian, can be any number when using gaussian

% Choose parameters to vary
% NO VARIATIONS
% The islands are all identical

varparameter = 0;
sigma_par = 0; % normalised standard deviation
mean_par = 1;
s1 = 1; s2 = 1; % upper and lower limit of integration same

% SHAPE
% The island volume is fixed as the island in-plane ellipticity (b/a) is
% varied. For other ways of varying shape, relevant sections in pswitch and
% demagdatagen and fielddatagen should be updated

% varparameter = 1;
% sigma_par = 0.075; % normalised standard deviation, i.e sigma_shape/mean shape
% mean_par = 1; % normalised mean shape (ellipticity)
% n1 = 10;    % number determining upper limit of integration, see x1
% n2 = 10;    % number determining lower limit of integration, see x2
% s1 = mean_par - n1*sigma_par;  % upper limit of intergration
% s2 = mean_par + n2*sigma_par;  % lower limit of intergration

% SIZE 
% The island in-plane ellipticity (b/a=1) and thickness is fixed as
% the volume is varied. For other ways of varying size, relevant sections
% in pswitch and demagdatagen and fielddatagen should be updated

% varparameter = 2;
% sigma_par = 0.075; % normalised standard deviation, i.e sigma_size/mean size
% mean_par = 1; % normalised mean size (volume)
% n1 = 10;    % number determining upper limit of integration, see x1
% n2 = 10;    % number determining lower limit of integration, see x2
% s1 = mean_par - n1*sigma_par;  % upper limit of intergration
% s2 = mean_par + n2*sigma_par;  % lower limit of intergration

% POSITION:JITTER
% The island position is varied from the expected position

% varparameter = 3;
% sigma_par = 0.075; % normalised standard deviation, i.e sigma_position/downtrackperiod
% mean_par = 1; %0 % normalised mean position
% n1 = 10;    % number determining upper limit of integration, see x1
% n2 = 10;    % number determining lower limit of integration, see x2
% s1 = mean_par - n1*sigma_par;  % upper limit of intergration
% s2 = mean_par + n2*sigma_par;  % lower limit of intergration


%LEAVE COMMENTED UNLESS VARYING JITTER TO CROSSTRACK 
% jitter_down_or_cross = 1; %position variatians are crosstrack

% CRYSTALLINE ANISOTROPY
% The crystalline anisotropy constant, k1, is varied from the expected
% value

% varparameter = 4;
% sigma_par = 0.075; % normalised standard deviation, i.e sigma_k1/mean k1
% mean_par = 1; % normalised mean k1
% n1 = 10;    % number determining upper limit of integration, see x1
% n2 = 10;    % number determining lower limit of integration, see x2
% s1 = mean_par - n1*sigma_par;  % upper limit of intergration
% s2 = mean_par + n2*sigma_par;  % lower limit of intergration

var_prop = [s1 s2 mean_par sigma_par varparameter jitter_down_or_cross w dist_type];


%% TOLERANCES FOR INTEGRATION

var_tol = [1e-8 1e-10]; % relative and abolute tolerance for integrating variation parameters: default [1e-6 1e-10] to disable set to zero
pswitch_rel_tol = 1e-6; % switching probability relative tolerance: default 1e-6
pswitch_abs_tol = 1e-6; % switching probability absolute tolerance: default 1e-10: to disable set to 0
headfield_tol = 1e-6;    % used in averaging head field over island volume
step_vect = [0 0; 0 0; 0 0]; % range of subset of the head field matrix for interpolation, all zeros means use all data
pswitch_quadrature_option = 1; % choose 0 for quad, 1 for quadgk
var_quadrature_option = 0; % choose 0 for quad, 1 for quadgk


tol_prop = [pswitch_rel_tol demag_tol headfield_tol pswitch_abs_tol pswitch_quadrature_option var_quadrature_option var_tol step_vect];

%% INTERPOLATION OF VARIATION PARAMETERS

% Initialising
interp_demag = 0; % for no interpolation in demagnetisation factors
interp_field = 0; % for no interpolation of volume averaged head fields
demag_data = 0; % for no interpolation of demagnetisation factors
h_data = 0; % for no interpolation of volume averaged head fields
x_data = 1; % for no interpolation of volume averaged head fields
y_data = 1; % for no interpolation of volume averaged head fields
s_data = 1; % for no interpolation of volume averaged head fields


% DEMAGNETISAION FIELD INTERPOLATION
% Note:Interpolation is NOT necessary for the prism case since expression is analytic
% 
npts =50; % number of points for interpolation
demagdatagen(a, b, t_h, alpha, islandgeo, varparameter, s1, s2, npts, demag_tol); % considering the hard layer only
demag_data = load('demagdata.m'); 
interp_demag = 1;


% HEAD FIELD INTERPOLATION

interp_field = 1; 
interp_prop = [interp_demag interp_field];

%% FILL VARIABLE ARRAY AND PRINT TO FILE 

t = now;
dir_name =strcat('Switching_Probability_Variables');
mkdir(dir_name); %create the directory
cd(dir_name);


filename = strcat('islandgeo_prop.m'); 
fid=fopen(filename,'a+');
fprintf(fid,'%g %g %g %g %g %g %g %g %g %g %g %g %g\n',[islandgeo a b t_h alpha islandposition_d islandposition_a v_h downtrackperiod  crosstrackperiod v_s area t_s]);

fclose(fid);

filename = strcat('islandmag_prop.m'); 
fid=fopen(filename,'a+');
fprintf(fid,'%g %g %g %g %g %g %g %g %g %g %g %g %g %g\n',[muo m_h h_k_h nxx_h nyy_h nzz_h nxx_s nyy_s nzz_s h_k_s, m_s, j_xc, s1_factor, s2_factor]);
fclose(fid);

filename = strcat('head_prop.m'); 
fid=fopen(filename,'a+');
fprintf(fid,'%g %g %g %g %g %g %g %g %g %g %g %g %g %g\n',[headtype hg phih gapsize polesize flyheight headposition_d headposition_a vel tau realheadposition_d realheadposition_a interlayer downtrack_travel]);
fclose(fid); 

filename = strcat('thermal_prop.m'); 
fid=fopen(filename,'a+');
fprintf(fid,'%g %g %g\n',[temp kb attfreq*write_attempts]);
fclose(fid); 

filename = strcat('tptiny_prop.m'); 
fid=fopen(filename,'a+');
fprintf(fid,'%g %g %g %g %g\n',[twait tptinysteps t_upper_sub tptiny_sub_steps tperiod]);
fclose(fid);  

filename = strcat('var_prop.m'); 
fid=fopen(filename,'a+');
fprintf(fid,'%g %g %g %g %g %g %g %g\n',[s1 s2 mean_par sigma_par varparameter jitter_down_or_cross w dist_type]);
fclose(fid);  

filename = strcat('tol_prop.m'); 
fid=fopen(filename,'a+');
fprintf(fid,'%g %g %g %g %g %g %g %g %g\n',[pswitch_rel_tol demag_tol headfield_tol pswitch_abs_tol pswitch_quadrature_option var_quadrature_option var_tol step_vect]);
fclose(fid);  

filename = strcat('interp_prop.m'); 
fid=fopen(filename,'a+');
fprintf(fid,'%g %g\n',[interp_demag interp_field]);
fclose(fid);  
%%
