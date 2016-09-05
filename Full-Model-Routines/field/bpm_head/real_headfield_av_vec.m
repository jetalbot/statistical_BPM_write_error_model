function [hx_av hy_av hz_av] = real_headfield_av_vec(realhead_pos_prop, realhead_field_prop, head_prop, islandgeo_prop, tol_prop, step_vect)
hg = head_prop(2);
islandgeo = islandgeo_prop(1);
a = islandgeo_prop(2);
b = islandgeo_prop(3);
t = islandgeo_prop(4);

x_vect = realhead_pos_prop(1,1):realhead_pos_prop(1,2):realhead_pos_prop(1,3);
y_vect = realhead_pos_prop(2,1):realhead_pos_prop(2,2):realhead_pos_prop(2,3);
z_vect = realhead_pos_prop(3,1):realhead_pos_prop(3,2):realhead_pos_prop(3,3);

hx = realhead_field_prop(:,:,:,1);
hy = realhead_field_prop(:,:,:,2);
hz = realhead_field_prop(:,:,:,3);

% to be verified: y-coordinate is down track, x-coordinate is cross track
y_coord = head_prop(11) + islandgeo_prop(6) - head_prop(7);  % actual y-coordinate, since at t=0, tsw=0, only first two terms contribute
x_coord = head_prop(12) + islandgeo_prop(7) - head_prop(8);  % actual x-coordinate, since at t=0, tsw=0, only first two terms contribute
z_coord = head_prop(13); % interlayer spacing
r_vec_center = [x_coord y_coord z_coord];

tol = tol_prop(3);

switch islandgeo
    case 1
        %         disp('truncated elliptic cone')
        alpha = islandgeo_prop(5);
        [hx_av hy_av hz_av] = headfieldcone_av_vec(a, b, t, alpha, r_vec_center, x_vect, y_vect, z_vect, step_vect, tol, hx, hy, hz);
        hx_av = hg.*hx_av;  % head field varies with time according to hg, i.e from 0 to 1 (saturation)
        hy_av = hg.*hy_av;
        hz_av = hg.*hz_av;

    case 2
        %         disp('elliptic cylinder')
        [hx_av hy_av hz_av] = headfieldcyl_av_vec(a, b, t, r_vec_center, x_vect, y_vect, z_vect, step_vect, tol, hx, hy, hz);
        hx_av = hg.*hx_av;  % head field varies with time according to hg, i.e from 0 to 1 (saturation)
        hy_av = hg.*hy_av;
        hz_av = hg.*hz_av;
    case 3
        %         disp('prism')
        [hx_av hy_av hz_av] = headfieldpris_av_vec(a, b, t, r_vec_center, x_vect, y_vect, z_vect, step_vect, tol, hx, hy, hz);
        hx_av = hg.*hx_av;  % head field varies with time according to hg, i.e from 0 to 1 (saturation)
        hy_av = hg.*hy_av;
        hz_av = hg.*hz_av;

    otherwise
        disp('Unknown geometry.')
end
end