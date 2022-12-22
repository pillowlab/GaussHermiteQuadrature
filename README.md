# GaussHermiteQuadrature
Gauss-Hermite Quadrature (for approximating integrals w.r.t. a Gaussian density)

------- 

**Theory**

The basic idea in Gauss-Hermite Quadrature is that we can evaluate an integral of the form of the product of some function $$f(x)$$ and a `weighting function' $$p(x)$$ as an appropriately weighted sum of function evaluations at a specified set of points.  In other words:

$$\int f(x) p(x) dx \approx \sum_{i=1}^N w_i f(x_i),$$

where the weights $$w_i$$ and evaluation points $$x_i$$ come from the theory of orthogonal polynomials.

