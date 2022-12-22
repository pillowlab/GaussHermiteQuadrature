% Test out using Gauss-Hermite Quadrature for evaluating numerical integral


% ====  Make Hermite polynomials of different orders ==========

n = 5;  % set order for Gauss-Hermite polynomial (higher -> more accurate)
[rr,ww] = compGaussHermiteQuadCoeffs(n); % get points and weights 

n2 = 7;  % set order for Gauss-Hermite polynomial (higher -> more accurate)
[rr2,ww2] = compGaussHermiteQuadCoeffs(n2); % get points and weights 

n3 = 10;  % set order for Gauss-Hermite polynomial (higher -> more accurate)
[rr3,ww3] = compGaussHermiteQuadCoeffs(n3); % get points and weights 

%%  ===== Set function to integrate ======

% Function we wish to integrate against a Gaussian density
fptr = @(x)(1./(1 + exp(-0.7*x-0.33)))*.8 + .1*sin(1.1*x);

% Gaussian to integrate over
mu = 5;  % mean of Gaussian
sigma = 3; % stdev of Gaussian

%% ===== Evaluate integral using Gauss-Hermite quadrature ========

fvals= fptr(rr*sigma + mu);  % evaluate function at these points
Fintegral = fvals'*ww; % evaluate integral using G-H quadrature

fvals2= fptr(rr2*sigma + mu);  % evaluate function at these points
Fintegral2 = fvals2'*ww2; % evaluate integral using G-H quadrature

fvals3= fptr(rr3*sigma + mu);  % evaluate function at these points
Fintegral3 = fvals3'*ww3; % evaluate integral using G-H quadrature


%% ====  Compute integral numerically using a grid ================

xrnge = mu + [-1 1]*sigma*10; % set range for numerical integral
nx = 1000; % number of grid points to use
dx = diff(xrnge)/nx; % grid spacing
xgrid = xrnge(1)+dx/2:dx:xrnge(2); % grid of points for evaluating func

px = normpdf(xgrid,mu,sigma);  % Gaussian density on grid
fx = fptr(xgrid);  % evaluate function on grid

% Evaluate the function numerically
Fnumerical = sum(fx.*px)*dx;

%%  Report results

fprintf('---------------------------------------------\n');
fprintf('Comparing Gauss-Hermite and Reimann integrals\n');
fprintf('---------------------------------------------\n');
fprintf('Reimann integral (%d points):      %.4f\n', nx, Fnumerical);
fprintf('Gauss-Hermite quadrature (order=%d):  %.4f (err=%8.4f)\n', n,Fintegral,Fnumerical-Fintegral);
fprintf('Gauss-Hermite quadrature (order=%d):  %.4f (err=%8.4f)\n', n2,Fintegral2,Fnumerical-Fintegral2);
fprintf('Gauss-Hermite quadrature (order=%d): %.4f (err=%8.4f)\n', n3,Fintegral3,Fnumerical-Fintegral3);

% Make plot showing function f(x) 
subplot(211);
plot(xgrid,px, xgrid,fx,rr3*sigma+mu, fvals3,'o');
set(gca,'xlim',xrnge);
legend('Gaussian density','f(x)','evaluation pts', 'location', 'northwest');
title('function f(x)');  
xlabel('x'); box off;

% Make plot showing N(mu, sig) and f(x)*N(mu,sig)
subplot(212); 
plot(xgrid,px, xgrid,fx.*px, rr3*sigma+mu, fvals3.*normpdf(rr3*sigma+mu, mu, sigma),'o');
set(gca,'xlim',xrnge);
legend('Gaussian density', 'f(x) * Gaussian', 'evaluation pts', 'location', 'northwest');
title('Gaussian and f(x) * Gaussian');
xlabel('x'); box off;