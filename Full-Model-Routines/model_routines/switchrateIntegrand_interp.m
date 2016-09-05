% Interpolate the sub steps, tp_dat, to find a probability, tprob_dat, within the
% switching probability array tp. This can then be integrated to find the
% switching proability. 
function tprob = switchrateIntegrand_interp(tp, tp_dat, tprob_dat)
tprob = interp1(tp_dat, tprob_dat, tp, 'pchip','extrap');
end