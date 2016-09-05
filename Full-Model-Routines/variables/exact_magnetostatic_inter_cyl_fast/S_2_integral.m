function out = S_2_integral(xo, a_1, a_2, m_upper)
summation =0;
m_upper
for m=0:m_upper
    m
    xo
    hyp = hypergeom([-m,-1-m],2,(a_2./a_1).^2).*((-(a_1./(2*xo)).^2).^m)   
    summation = summation + gamma(1 + 2*m)./(factorial(m).*gamma(m+2)).*hypergeom([-m,-1-m],2,(a_2./a_1).^2).*((-(a_1./(2*xo)).^2).^m);
end
out = a_1*a_2*(1/4)*xo.^(-1).*summation;
end