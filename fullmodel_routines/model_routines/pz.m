%% TRUNCATED GAUSSIAN PIECEWISE FUNCTION
% Calculates a piecewise Gaussian function up to truncated tails. The shape
% and positon of these tails are determined by $n*sigma$, where $n$ is a
% multiple value chosen so that the area of interest is in, in this case,
% the tail of the distribution, and $sigma$ is the standard deviation. A
% point $x'$ is then chosen so that the area under the curve remains the
% same as the Gaussian. 
%%
function out = pz(s,var_prop,n)

% Get initial values
mean = var_prop(1); 
sigma = var_prop(2);

% The point at which we deviate from the Gaussian, to the left (negative
% side of the mean) and the right (positive side of the mean).
xleft = mean - n*sigma; % starting point of truncation for the negative side
xright = mean + n*sigma; % starting point of truncation for the positive side


% Find the x' value, where the truncation ends. 
% Get the integral of the Gaussian between 2*sigma and -2*sigma
k = @(x)gauss(x,var_prop);
gaussnorm_2s = quadgk(k,xleft, xright);
% Get the coordinates for that starting point of the truncation and hence
% find the coordinate for x' that will make the area of the entire
% distribution the same as a Gaussian distribution.
y1 = normpdf(n*sigma+mean,mean,sigma);
x1 = (1-gaussnorm_2s)/y1;

% x' coordinates
xneg = -x1+xleft;
xpos = x1+xright;
             
% Get the equations for the truncating line                
M1 = y1/x1; %gradient, negative side
M2 = -M1; %gradient, positive side
C1 = - M1*xneg; %interception
C2 = - M2*xpos; %interception

%Piecewise function of the truncated distribution. For comparision the
%Gaussian distribution in the same form is 
% pz1 =(x < xneg).* normpdf(x,mean,sigma);
% pz2 =((x > xneg)&(x<= xleft)).*(normpdf(x,mean,sigma));
% pz3 = ((x > xleft)&(x<= xright)).*(normpdf(x,mean,sigma));
% pz4 = ((x > xright)&(x<= xpos)).*(normpdf(x,mean,sigma));
% pz5 = (x > xpos).*(normpdf(x,mean,sigma));
% out = pz1+pz2+pz3+pz4+pz5;

pz1 =((s < xneg).*(0));
pz2 =((s > xneg)&(s<= xleft)).*((M1*s)+C1);
pz3 = ((s > xleft)&(s<= xright)).*(normpdf(s,mean,sigma));
pz4 = ((s > xright)&(s<= xpos)).*((M2*s)+C2);
pz5 = (s > xpos).*(0);
out = pz1+pz2+pz3+pz4+pz5;
