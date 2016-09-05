function [havx havy havz] = karlqfield_av(head_prop, islandgeo_prop, tol_prop)
hg = head_prop(2);
phih = head_prop(3);
gapsize = head_prop(4);
polesize = head_prop(5);
zh = head_prop(6);
rh = head_prop(7)- islandgeo_prop(6);  % only down track separation in this case

a = islandgeo_prop(2);
b = islandgeo_prop(3);
t = islandgeo_prop(4);

tol = tol_prop(3);
switch islandgeo_prop(1)
    case 1
%         disp('truncated elliptic cone')
        alpha = island_prop(5);
        [havx havy havz] = karlqfieldcone_av(hg, phih, gapsize, polesize, rh, zh, a, b, t, alpha, tol);
    case 2
%         disp('elliptic cylinder')
        [havx, havy, havz] = karlqfieldcyl_av(hg, phih, gapsize, polesize, rh, zh, a, b, t, tol);
    case 3
%         disp('prism')
        d = gapsize/2 -zh -t;
        [havx havy havz] = karlqfieldpris_av(hg, phih, gapsize, polesize, d, rh, a, b, t, tol);
    otherwise
        disp('Unknown variation parameter.')
end

end