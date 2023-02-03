function [rr,ww] = compGaussHermiteQuadCoeffs0(n)
% [rr,ww] = compGaussHermiteQuadCoeffs0(n)
% 
% Compute roots and weights needed for Gauss-Hermite Quadrature using a
% standard normal weighting function. 
%
% This function uses a recursive formula to get coefficients; the function 
% 'compGaussHermiteQuadCoeffs.m' instead computes them using closed-form
% expressions (and should thus be slightly more efficient).   
%
% Approximates definite integral: 
%
%      F  = \int_-infy^infty f(x) N(x) dx 
%  as
%
%      F_approx = \sum_{i=1}^n ww_i f(rr_i)
%
% where N(x) = 1/sqrt(2pi) exp(-x^2/2) and f(x) is any smooth function 
%
%
% INPUT
% -----
% n [1 x 1] - order of polynomial (integer >= 1). (higher -> more accurate)
%
% OUTPUT
% ------
% rr [n x 1] - roots (location to evaluate function)
% ww [n x 1] - weights for Gauss-Hermite Quadrature
%
%
% ---------------
% Note: this differs slightly from standard Gauss-Hermite Quadrature, in
% which the weighting function N(x) is exp(-x^2)
%
% To approximation to integrals w/ non-standard normal densities, i.e.,   
%    F = \int f(x) N(x ; mu, sigma^2)
%  use   
%    F_approx = sum ww_i f(rr_i * sigma + mu)


% test if n is an integer >= 1
if abs(mod(n,1))>1e-10 || n<1
    error('polynomial order ''n'' must be an integer >= 1');
end

% allocate space for polynomial coeffs
cc = zeros(n+1,n+1); 

% initialize first two rows
cc(1:2,1)=[1;2]; 

% Construct polynomial coeffs recursively
for k=2:n
   cc(k+1,1:k+1)=2*[cc(k,1:k) 0]-2*(k-1)*[0 0 cc(k-1,1:k-1)];
end
Hmcoefs = cc(n,1:n);  % coeffs for order n-1 polynomial
Hncoefs = cc(n+1,:);  % coeffs for order n polynomial

% Find roots (points at which to evaluate polynomial)
Hzeros=roots(Hncoefs);

% Compute Coeficients
num = (2.^(n-1)*(factorial(n-1)))/n; % numerator
pvals = polyval(Hmcoefs,Hzeros).^2;  % denominator
ww = num./pvals;  % weights

% transform roots by sqrt(2) to obtain points at which to evaluate func
rr = Hzeros*sqrt(2);
