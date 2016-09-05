%% GAUSSIAN PIECEWISE FUNCTION
% Calculates a piecewise Gaussian function. 
%% 
function out = px(s,var_prop,n)

% Get initial values
mean = var_prop(1); 
sigma = var_prop(2);

% The point at which we deviate from the Gaussian, to the left (negative
% side of the mean) and the right (positive side of the mean).
xleft = mean - n*sigma; % starting point of truncation for the negative side
xright = mean + n*sigma; % starting point of truncation for the positive side


% Find the x' value, where the truncation ends. 
% Get the integral of the Gaussian between 2*sigma and -2*sigma
k = @(x)normpdf(x,mean,sigma);
gaussnorm_2s = quadgk(k,xleft, xright);
% Get the coordinates for that starting point of the truncation and hence
% find the coordinate for x' that will make the area of the entire
% distribution the same as a Gaussian distribution.
y1 = normpdf(n*sigma+mean,mean,sigma);
x1 = (1-gaussnorm_2s)/y1;

% x' coordinates
xneg = -x1+xleft;
xpos = x1+xright;
             
%Piecewise function of the distribution.
pz1 =(s < xneg).* normpdf(s,mean,sigma);
pz2 =((s > xneg)&(s<= xleft)).*(normpdf(s,mean,sigma));
pz3 = ((s > xleft)&(s<= xright)).*(normpdf(s,mean,sigma));
pz4 = ((s > xright)&(s<= xpos)).*(normpdf(s,mean,sigma));
pz5 = (s > xpos).*(normpdf(s,mean,sigma));
out = pz1+pz2+pz3+pz4+pz5;
end