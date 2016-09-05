function q = myquad(f,a,b,tol,trace,varargin)
q = quadgk(@(x)f(x,varargin{:}),a,b,'AbsTol',tol);
end
