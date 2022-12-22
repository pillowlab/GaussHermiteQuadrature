# GaussHermiteQuadrature
Gauss-Hermite Quadrature (for approximating integrals w.r.t. a Gaussian density)

------- 

###Gauss-Hermite Quadrature###

The basic idea in Gauss-Hermite Quadrature is that we can evaluate an integral of the form of the product of some function $f(x)$ and a `weighting function' $p(x)$ as an appropriately weighted sum of function evaluations at a specified set of points.  In other words:

$$\int_{-\infty}^\infty f(x) p(x) dx \;\approx\; \sum_{i=1}^N w_i f(r_i),$$

where $p(x) = \exp(-x^2)$, and the weights $w_i$ and evaluation points $r_i$ come from the theory of Hermite polynomials (which are orthogonal polynomials w.r.t. the weighting function $p(x)$).  

In THIS repository, however, the weighting function used is a standard normal Gaussian density, $p(x) = 1/\sqrt(2\pi) \exp(x^2 / 2)$.  

------- 

### Using this repository ###

-  Call *compGaussHermiteQuadCoeffs* to obtain the evaluation points $r_i$ and the weights $w_i$

> [rr,ww] = compGaussHermiteQuadCoeffs(n);


