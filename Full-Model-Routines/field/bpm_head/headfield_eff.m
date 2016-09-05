function h_eff = headfield_eff(hx, hy, hz)
h = sqrt(hx.^2 + hy.^2 + hz.^2);
thetah = acos(hz./h);
h_eff = h./((cos(thetah).^2).^(1/3)+ (sin(thetah).^2).^(1/3)).^(-3/2);
end