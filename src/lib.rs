/// bsfun - basis spline fun(ctions)
/// 2023-05-11 Dorian Prill
/// A very simple library for evaluating NURBS (non-uniform rational b-spline) 
/// curves and surfaces in arbitrary dimension.
/// Consists of only three functions, no structures or traits.:
///     bspline_basis()
///     rational_bspline_basis()
///     nurbs_curve_point()


// Calculates the B-Spline Basis Function of `degree` at index `i` 
// over the non-uniform `knots`, evaluated at `t`.
// A degree zero B-Spline is defined as
//  Bi,0(t) :=  1 if ti <= t < ti+1
//              0 otherwise
//  allowing for the recursive definition of a degree k B-Spline as
//  Bi,k(t) := (t - ti) / (ti+k - ti) * Bi,k-1(t) + (ti+k+1 - t) / (ti+k+1 - ti+1) * Bi+1,k-1(t)
pub fn bspline_basis(i: usize, degree: usize, knots: &[f64], t: f64) -> f64 {
    if degree == 0 {
        if (t >= knots[i] && t < knots[i+1]) || (i+1 == knots.len()-1 && t == knots[i+1]) { 1.0 } 
        else { 0.0 }
    } else {
        let a = if (knots[i+degree] - knots[i]) == 0.0 { 0.0 } else { (t - knots[i]) / (knots[i+degree] - knots[i]) * bspline_basis(i, degree-1, knots, t) };
        let b = if (knots[i+degree+1] - knots[i+1]) == 0.0 { 0.0 } else { (knots[i+degree+1] - t) / (knots[i+degree+1] - knots[i+1]) * bspline_basis(i+1, degree-1, knots, t) };
        a + b
    }
}

// Calculates the rational (weighted) B-Spline Basis Function of `degree` 
// at index `i` over the `knots`, evaluated at `t`.
pub fn rational_bspline_basis(i: usize, degree: usize, knots: &[f64], weights: &[f64], t: f64) -> f64 {
    let numerator = weights[i] * bspline_basis(i, degree, knots, t);
    let denominator: f64 = (0..weights.len()).map(|j| weights[j] * bspline_basis(j, degree, knots, t)).sum();
    if denominator.abs() < std::f64::EPSILON { 0.0 } else { numerator / denominator }
}

// Calculates a point on a NURBS curve
// This function assumes that the weights slice and the control_points slice have 
// the same length, and that all control points have the same dimension. 
// It also assumes that the knots slice has the correct length for the number 
// of control points and the degree of the curve. 
// If these assumptions aren't met, the function may panic or return incorrect results.
pub fn nurbs_curve_point<'a>(t: f64, control_points: &'a [Vec<f64>], weights: &'a [f64], knots: &'a [f64], degree: usize) -> Vec<f64> {
    let n_dims = control_points[0].len();
    let mut point = vec![0.0; n_dims];

    for i in 0..control_points.len() {
        let basis = rational_bspline_basis(i, degree, knots, weights, t);
        for k in 0..n_dims {
            point[k] += basis * control_points[i][k];
        }
    }

    point
}



// Calculates a point on a NURBS surface
// Same assumptions as for nurbs_curve_point()
pub fn nurbs_surface_point<'a>(u: f64, v: f64, control_points: &'a [Vec<Vec<f64>>], weights: &'a [Vec<f64>], knots_u: &'a [f64], knots_v: &'a [f64], degree_u: usize, degree_v: usize) -> Vec<f64> {
    let n_dims = control_points[0][0].len();
    let mut point = vec![0.0; n_dims];

    for i in 0..control_points.len() {
        for j in 0..control_points[i].len() {
            let basis_u = rational_bspline_basis(i, degree_u, knots_u, &weights[i], u);
            let basis_v = rational_bspline_basis(j, degree_v, knots_v, &weights[j], v);
            let basis = basis_u * basis_v;

            for k in 0..n_dims {
                point[k] += basis * control_points[i][j][k];
            }
        }
    }

    point
}





#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bspline_basis_degree_zero() {
        let knots = vec![0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0];
        for i in 0..knots.len()-1 {
            for t in 0..=100 {
                let t = t as f64 / 100.0;
                let result = bspline_basis(i, 0, &knots, t);
                dbg!(i, t, result);
                if t >= knots[i] && t < knots[i+1] {
                    assert!((result-1.0).abs() < std::f64::EPSILON);
                } else {
                    assert!(result.abs() < std::f64::EPSILON);
                }
            }
        }
        
    }

    #[test]
    fn test_bspline_basis_degree_one() {
        // simple test by checking that they are 1 at the midpoint of 
        // their knot span and 0 outside of their knot span, 
        //since we know that they linearly interpolate between these values
        let knots = vec![0.0, 0.0, 1.0, 2.0, 3.0, 3.0];

        for i in 0..knots.len()-2 {
            for t in 0..=100 {
                let t = t as f64 / 100.0;
                let result = bspline_basis(i, 1, &knots, t);
                dbg!(i, t, result);
                if t >= knots[i] && t < knots[i+1] {
                    assert!(result > (1.0 - (knots[i+1] - t) - 0.000001) && result < (1.0 - (knots[i+1] - t) + 0.000001));
                } else if t >= knots[i+1] && t < knots[i+2] {
                    assert!(result > ((knots[i+2] - t) - 0.000001) && result < ((knots[i+2] - t) + 0.000001));
                } else {
                    assert!(result > (-0.000001) && result < 0.000001);
                }
            }
        }
    }


    #[test]
    fn test_rational_bspline_basis_degree_zero() {
        let knots = vec![0.0, 0.0, 1.0, 2.0, 3.0, 3.0];
        let weights = vec![1.0, 1.0, 1.0, 1.0, 1.0];

        for i in 0..knots.len()-1 {
            for t in 0..=100 {
                let t = t as f64 / 100.0;
                let result = rational_bspline_basis(i, 0, &knots, &weights, t);
                dbg!(i, t, result);
                if t >= knots[i] && t < knots[i+1] {
                    assert!(result > (1.0-0.000001) && result < (1.0+0.000001));
                } else {
                    assert!(result > (-0.000001) && result < 0.000001);
                }
            }
        }
    }


    #[test]
    fn test_rational_bspline_basis_degree_one() {
        // This test function computes the total weighted value of all the basis functions, 
        // and uses it to compute the expected value of each rational B-spline basis function
        let knots = vec![0.0, 0.0, 1.0, 2.0, 3.0, 3.0];
        let weights = vec![1.0, 2.0, 3.0, 4.0];
        // the loop goes up to knots.len() - 3 instead of knots.len() - 2. 
        // This is because for degree 1 basis functions, we need to 
        // consider one extra knot beyond the current knot span 
        // when computing the basis function value.
        for i in 0..knots.len()-2-1 {
            for t in 0..=100 {
                let t = t as f64 / 100.0;
                let result = rational_bspline_basis(i, 1, &knots, &weights, t);
                dbg!(i, t, result);
                let mut total = 0.0;
                for j in 0..weights.len() {
                    total += weights[j] * bspline_basis(j, 1, &knots, t);
                }
                if total.abs() > 0.000001 {
                    let expected = (weights[i] * bspline_basis(i, 1, &knots, t)) / total;
                    assert!(result > (expected - 0.000001) && result < (expected + 0.000001));
                } else {
                    assert!(result > (-0.000001) && result < 0.000001);
                }
            }
        }
    }


}

