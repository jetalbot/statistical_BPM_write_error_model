function h  = sField(theta)
t = nthroot(tan(theta), 3);
a = (1 - t.^2 + t.^4).^(1/2);
b = 1 + t.^2;
h = real(a./b); % real and positive, only dealing with magnitudes
end
