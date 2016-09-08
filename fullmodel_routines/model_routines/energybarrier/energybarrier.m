function [barrier1 barrier2] = energybarrier(hv, phiv)
% This functions returns the barrier ideally for h positive, to find the other barrier replace h by
% -h
phivec =ones(size(hv)).*phiv;
barrier1 = zeros(size(hv));
barrier2 = zeros(size(hv));

for j=1:length(hv)
    h=hv(j);
    phi=phivec(j);
    if isnan(h) || isnan(phi) % this occurs hx=hy=hz=0, polar and azimuthal angles become indeterminate => heff becomes indeterminate too
        barrier1(j) = 0.5;
        barrier2(j) = 0.5;
    elseif sin(phi)==0 || sign(h)==0 % here the maxima and minima become degenerate and simpler to do it this way
        if abs(h) <= sField(phi) || sign(h)==0
            barrier1(j) = ((1 + h).^2)/2;
            barrier2(j) = ((1 - h).^2)/2;
        else
            barrier1(j) = 2.*abs(h);
            barrier2(j) = 2.*abs(h);
        end
    else
        v = (h.^2.*sin(phi).*cos(phi) + sqrt((h.^2.*sin(phi).*cos(phi)).^2 + ((h.^2 -1)/3).^3 )).^2;
        w = (h.^2.*sin(phi).*cos(phi) - sqrt((h.^2.*sin(phi).*cos(phi)).^2 + ((h.^2 -1)/3).^3 )).^2;

        vcrt = nthroot(real(v), 3).*(real(v)==v) + v.^(1/3).*(real(v)~=v);
        wcrt = nthroot(real(w), 3).*(real(w)==w) + w.^(1/3).*(real(w)~=w);

        t = vcrt + wcrt;

        c1 = sqrt(t + h.^2.*sin(phi).^2 - 2*(h.^2 -1)/3);
        den1 = sqrt(t + h.^2.*sin(phi).^2 - 2*(h.^2 -1)/3);
        den2 = sqrt(t + h.^2.*cos(phi).^2 - 2*(h.^2 -1)/3);

        den1 = eps*(abs(den1)<=eps) + den1;
        den2 = eps*(abs(den2)<=eps) + den2;

        srt1 = sqrt(2*h.^2.*sin(phi).^2 - 4*(h.^2 -1)/3 - t + 2*h.*sin(phi).*(h.^2.*cos(phi).^2 + 1)./den1);
        srt2 = sqrt(2*h.^2.*cos(phi).^2 - 4*(h.^2 -1)/3 - t + 2*h.*cos(phi).*(h.^2.*sin(phi).^2 + 1)./den2);

        barrier1(j) = (c1 + h.*sin(phi)).*srt1/2 + h.*cos(phi).*srt2;

        h =-h;

        srt1 = sqrt(2*h.^2.*sin(phi).^2 - 4*(h.^2 -1)/3 - t + 2*h.*sin(phi).*(h.^2.*cos(phi).^2 + 1)./den1);
        srt2 = sqrt(2*h.^2.*cos(phi).^2 - 4*(h.^2 -1)/3 - t + 2*h.*cos(phi).*(h.^2.*sin(phi).^2 + 1)./den2);

        barrier2(j) = (c1 + h.*sin(phi)).*srt1/2 + h.*cos(phi).*srt2;
    end
end

end
