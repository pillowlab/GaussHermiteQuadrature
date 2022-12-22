# GaussHermiteQuadrature
Gauss-Hermite Quadrature (for approximating integrals w.r.t. a Gaussian density)

------- 

### Gauss-Hermite Quadrature ###

The basic idea in Gauss-Hermite Quadrature is that we can evaluate an integral of the form of the product of some function $f(x)$ and a `weighting function' $p(x)$ as an appropriately weighted sum of function evaluations at a specified set of points.  In other words:

$$\int_{-\infty}^\infty f(x) p(x) dx \quad \approx \quad \sum_{i=1}^N w_i f(r_i),$$

where $p(x) = \exp(-x^2)$, and the weights $w_i$ and evaluation points $r_i$ come from the theory of Hermite polynomials (which are orthogonal polynomials w.r.t. the weighting function $p(x)$).  

In _THIS_ repository, however, the weighting function used is a standard normal Gaussian density, $p(x) = 1/\sqrt{2\pi} \exp(x^2 / 2)$, which means that the sum will evaluate $f(x)$ times a normal density instead of $\exp(-x^2)$.

------- 

### Using this repository ###

-  Set the order of the Hermite polynomial $n$.  (This must be an integer >= 1;  Higher $n$ will give higher accuracy, though requires more function evalauations).  

> n = 10;         % order for Gauss-Hermite polynomial

-  Call *compGaussHermiteQuadCoeffs* to obtain the evaluation points $r_i$ and the weights $w_i$

> [rr,ww] = compGaussHermiteQuadCoeffs(n);    % get points and weights
> 

-  Then, for any function of interest $f$, evaluate the integral of $f$ times a standard normal pdf over the reals via: 

> fIntegral = sum( f(rr) .* ww);  

- To evaluate the integral of $f$ times a Gaussian with mean $\mu$ and variance $\sigma^2$, instead use:

> fIntegral = sum( f(rr * sig + mu) .* ww));


- See the script *testGaussHermiteQuadrature* for an example application