% Test out using Gauss-Hermite Quadrature for evaluating numerical integral


% ====  Make Hermite polynomial ==========

% Compute roots and weights of polynomial

n = 25;  % order for Gauss-Hermite polynomial (higher -> more accurate)

[rr,ww] = compGaussHermiteQuadCoeffs(n);
toc;

tic;
[rr0,ww0] = compGaussHermiteQuadCoeffs0(n);
toc;

%%  ===== Set function to integrate ======

% parameters of ideal observer model
obs_b = 1.1;  % offset
obs_sig = 1; % internal noise stdev
q0 = 1; % reward for left choice
q1 = 2; % reward for right choice

% Function pointer for observer policy (i.e., function to integrate)
fptr = @(x)(1./(1 + exp(q0 - (q0+q1)*normcdf((x-obs_b)/obs_sig))));

% Gaussian weighting function to integrate over
stim = 1.5;  % stimulus (point on the x axis where we wish to evaluate psychometric function)
mu = stim;   % mean of Gaussian
sigma = obs_sig; % stdev

%% ===== Evaluate integral using Gauss-Hermite quadrature ========

fvals= fptr(rr*sigma + mu);  % evaluate function at these points
Fgausshermite = fvals'*ww; % evaluate integral using G-H quadrature

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
fprintf('Gauss-Hermite quadrature (order=%d): %.4f\n', n,Fgausshermite);
fprintf('Reimann integral:                    %.4f\n', Fnumerical);
fprintf('error = %f\n', Fnumerical-Fgausshermite);

% Make plot
plot(xgrid,fx,xgrid,px, xgrid,fx.*px);
set(gca,'xlim',xrnge);
legend('f(x)', 'Gaussian density', 'Gaussian * f(x)', 'location', 'northwest');
