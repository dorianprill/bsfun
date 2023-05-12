# bsfun

## Basis Spline Fun(ctions)

This is a super simple Rust library for working with basis splines and [NURBS (Non-Uniform Rational B-Splines)](https://en.wikipedia.org/wiki/Non-uniform_rational_B-spline) in arbitrary dimension with zero dependecies (besides `alloc::Vec`).

These five basic functions are used in `NURBSCurve` and `NURBSSurface` evaluation:

1. `bspline_basis()` Calculate the value of a basis spline at a given `t`. This is the basic building block for all splines and the rational basis functions.
2. `rational_bspline_basis_curve()` Calculate the value of a rational basis spline for curves at a given `t`
3. `rational_bspline_basis_surface()` Calculate the value of a rational basis spline for surfaces at a given value pair `(u, v)`
4. `nurbs_curve_point()` Calculate the value of a NURBS curve at a given value `t`
5. `nurbs_surface_point()` Calculate the value of a NURBS surface at a given value pair `(u, v)`

You could use `1`, `2` and `3` to build your own Splines. `4` and `5` are just functions for calculating the value of a NURBS curve or surface for a set of parameters directly.

> Caution: `4` and `5` are mostly for playing around and will not validate parameters for you. You have to make sure that the parameters are in the correct range and of correct shape yourself. If you don't want that use the structs `NURBSCurve` and `NURBSSurface` and their `eval()` methods

Since this library does not have enhanced functionality beyond curve evaluations, you may want to use a more practical NURBS implementation from e.g. the library `truck`.

## Visualization of B-Spline Basis Functions up to Degree 4

made with [plotters](https://github.com/plotters-rs/plotters)

![B-Spline Basis Degree 1](basis_degree_1.png)

![B-Spline Basis Degree 2](basis_degree_2.png)

![B-Spline Basis Degree 3](basis_degree_3.png)

![B-Spline Basis Degree 4](basis_degree_4.png)
