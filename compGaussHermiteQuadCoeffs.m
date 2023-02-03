function [rr,ww] = compGaussHermiteQuadCoeffs(n)
% [rr,ww] = compGaussHermiteQuadCoeffs(n)
% 
% Compute roots and weights needed for Gauss-Hermite Quadrature using a
% standard normal weighting function.  
%
% Approximates definite integral: 
%
%      F  = \int f(x) N(x) dx    (from -inf to inf)
%  as
%
%      F_approx = \sum_{i=1}^n ww_i f(rr_i)
%
% where N(x) = 1/sqrt(2pi) exp(-x^2/2) and f(x) is any smooth function 
%
% INPUT
% -----
% n [1 x 1] - order of polynomial (integer >= 1). (higher -> more accurate)
%
% OUTPUT
% -----
% rr [n x 1] - roots (location to evaluate function)
% ww [n x 1] - weights for Gauss-Hermite Quadrature
%
% ----------
% Note: this code differs slightly from standard Gauss-Hermite Quadrature, 
% in which the weighting function N(x) is exp(-x^2). Here we use a standard
% normal density instead.
%
% To approximation to integrals w/ non-standard normal densities, i.e.,   
%    F = \int f(x) N(x ; mu, sigma^2)
%  use   
%    F_approx = sum ww_i f(rr_i * sigma + mu)

% test if n is an integer >= 1
if abs(mod(n,1))>1e-10 || n<1
    error('polynomial order ''n'' must be an integer >= 1');
end

% Compute the log coefs for n'th Hermite polynomial Hn
iin = 0:floor(n/2); % indices such that n-2*ii are nonzero
Hnsigns = zeros(1,n+1);  % signs of coefs
Hnsigns(n-(n:-2:0)+1) = (-1).^iin; % compute signs of coeffs
logHncoefs = zeros(1,n+1); % log-coeffs 
logHncoefs(n-(n:-2:0)+1) = gammaln(n+1) - gammaln(iin+1) - gammaln(n-2*iin+1) + (n-2*iin)*log(2);

% Find roots of Hn (points at which to evaluate function)
Hzeros = real(roots(Hnsigns.*exp(logHncoefs-max(logHncoefs))));

% Compute the log coefs for n-1'th Hermite polynomial Hm
m = n-1;  % n-1
iim = 0:floor(m/2); % indices such that n-2*ii are nonzero
Hmsigns = zeros(1,m+1);  % signs of coefs
Hmsigns(m-(m:-2:0)+1) = (-1).^iim; % compute signs of coefs
logHmcoefs = zeros(1,m+1); % log-coeffs 
logHmcoefs(m-(m:-2:0)+1) = gammaln(m+1) - gammaln(iim+1) - gammaln(m-2*iim+1) + (m-2*iim)*log(2);
maxlogHmcoefs = max(logHmcoefs,[],2);  % max of each row
Hmcoefs = Hmsigns.*exp(logHmcoefs-maxlogHmcoefs);

% Compute weights
lognum = (n-1)*log(2)+gammaln(n)-log(n); % numerator
logdenom = 2*log(abs(polyval(Hmcoefs,Hzeros)));  % denominator
ww = exp(lognum-logdenom-2*maxlogHmcoefs);  % make weights

% transform roots by sqrt(2) to obtain points at which to evaluate func
rr = Hzeros*sqrt(2);
