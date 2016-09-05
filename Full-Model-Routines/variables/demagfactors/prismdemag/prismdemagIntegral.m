function [nxx nyy nzz] = prismdemagIntegral(a, b, c)
nzz = prismdemag(a, b, c);
nyy = prismdemag(c, a, b);
nxx = prismdemag(b, c, a);

end