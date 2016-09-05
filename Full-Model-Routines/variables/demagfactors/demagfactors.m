function [ nxx nyy nzz ] = demagfactors(a, b, t, alpha, islandgeo, demag_tol)
switch islandgeo
    case 1
%         disp('truncated elliptic cone')
        [ nxx nyy nzz ] = conedemagIntegral(a, b, t, alpha, demag_tol);
    case 2
%         disp('elliptic cylinder')
        [ nxx nyy nzz ] = cyldemagIntegral(a, b, t, demag_tol);
    case 3
%         disp('prism')
        [nxx nyy nzz] = prismdemagIntegral(a, b, t/2); % since c =t/2
    otherwise
        disp('Unknown geometry.')
end
end