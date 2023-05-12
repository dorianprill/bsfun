# bsfun

## Basis Spline Fun(ctions)

This is a super simple Rust library for working with basis splines and [NURBS (Non-Uniform Rational B-Splines)](https://en.wikipedia.org/wiki/Non-uniform_rational_B-spline) with zero dependecies (besides `alloc::Vec`).

It consists of only four functions:

1. `bspline_basis()` Calculate the value of a basis spline at a given `t`
2. `rational_bspline_basis()` Calculate the value of a rational basis spline at a given `t`
3. `nurbs_curve_point()` Calculate the value of a NURBS curve at a given value `t`
4. `nurbs_surface_point()` Calculate the value of a NURBS surface at a given value pair `(u, v)`

You could use `1` and `2` to build your own Splines. `3` and `4` are just functions for calculating the value of a NURBS curve or surface for a set of parameters.

Since this library does not provide any structs or traits you have take care of the proper shapes of knots, weights and control points for any given degree yourself. You may find a more practical NURBS implementation in the libraries `truck` and `capstan`.

## Visualization of B-Spline Basis Functions up to Degree 4

made with [plotters](https://github.com/plotters-rs/plotters)

![B-Spline Basis Degree 1](basis_degree_1.png)

![B-Spline Basis Degree 2](basis_degree_2.png)

![B-Spline Basis Degree 3](basis_degree_3.png)

![B-Spline Basis Degree 4](basis_degree_4.png)
